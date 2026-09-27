# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""The hand loop on a group of ground-truth points: fit, normal from neighbours' positions, tilt, refit.

The group is either every finite point within ``--radius`` metres of
``--centre`` whose normal is within ``--max-normal-deg`` of the centre's, or an
explicit ``--points`` list (a comma-separated list, or ``@file``). Each track
starts at the ground-truth position with the ``--init`` normal: ``gt``,
``meanview`` (the mean direction to its cameras), or ``tiltNN`` (the ground
truth's tilted NN degrees in a random direction). Each round every point's
normal is taken from the fitted positions of the group points within
``--nb-halves`` half-sizes of it (a plane through them; where they only span a
line, only the tilt that line determines is changed), then each track is
tilted to it and fitted again. Reported per round, against the ground truth:
the normal error, the position error in half-sizes, and the median ZNCC, over
the group and over its sloped points (normal 5 degrees or more from vertical).

    pixi run -e test python scripts/track_at_pixel/surface/neighbour_positions.py DATASET --centre 112 --radius 3 --init meanview
"""

import argparse

import numpy as np

from _common import UP, ang, load, unit


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dataset")
    ap.add_argument("--centre", type=int)
    ap.add_argument("--radius", type=float, default=1.5)
    ap.add_argument("--max-normal-deg", type=float, default=90.0)
    ap.add_argument("--points", help="comma-separated point ids, or @file")
    ap.add_argument("--init", default="meanview", help="gt, meanview or tiltNN")
    ap.add_argument("--rounds", type=int, default=4)
    ap.add_argument("--nb-halves", type=float, default=2.0)
    args = ap.parse_args()

    from sfmtool._sfmtool import bench as B

    from candidates.common import median_zncc, patch_normal

    rng = np.random.default_rng(0)
    ds, edited = load(args.dataset)
    if args.points:
        text = args.points
        if text.startswith("@"):
            text = open(text[1:]).read()
        group = [int(x) for x in text.replace("\n", ",").split(",") if x.strip()]
    else:
        n0 = ds.point_normal[args.centre]
        group = [
            int(p)
            for p in np.flatnonzero(ds.point_w != 0)
            if np.linalg.norm(ds.point_xyz[p] - ds.point_xyz[args.centre]) < args.radius
            and ang(ds.point_normal[p], n0) < args.max_normal_deg
        ]
    X_gt = ds.point_xyz[group]
    half = ds.point_half[group]
    d = np.linalg.norm(X_gt[:, None] - X_gt[None], axis=2) + np.eye(len(group)) * 1e9
    sloped = np.asarray([ang(ds.point_normal[p], UP) >= 5 for p in group])
    print(
        f"{len(group)} points ({sloped.sum()} sloped); nearest-neighbour spacing / "
        f"patch diameter: median {np.median(d.min(1) / (2 * half)):.2f}; views per "
        f"point median {np.median([len(ds.point_images[p]) for p in group]):.0f}"
    )

    tracks, n_gt = [], []
    for p in group:
        _, t = B.create_track(B.Bench(), edited, p)
        cam = ds.cameras[int(t.observations[0]["image"])]
        g = unit(patch_normal(t))
        n_gt.append(g if g @ (cam.center - ds.point_xyz[p]) > 0 else -g)
        tracks.append(t)
    n_gt = np.asarray(n_gt)

    def to_cam(t):
        c = ds.cameras[int(t.observations[0]["image"])].center
        return unit(c - np.asarray(t.position))

    def init_normal(i, t):
        if args.init == "gt":
            return n_gt[i]
        if args.init == "meanview":
            X = np.asarray(t.position)
            return unit(
                sum(
                    unit(ds.cameras[int(o["image"])].center - X) for o in t.observations
                )
            )
        if args.init.startswith("tilt"):
            deg = float(args.init[4:])
            r = unit(np.cross(n_gt[i], rng.normal(size=3)))
            return unit(np.cos(np.radians(deg)) * n_gt[i] + np.sin(np.radians(deg)) * r)
        raise ValueError(args.init)

    def fit(t):
        try:
            t, _ = B.fit(t, edited, ds.pyramids)
        except ValueError:
            pass
        return t

    def tilt(t, n):
        n = n if n @ to_cam(t) > 0 else -n
        try:
            t, _ = B.tilt_patch(t, edited, tuple(float(x) for x in n))
        except ValueError:
            pass
        return t

    def report(label, ts):
        X = np.asarray([t.position for t in ts])
        ne = np.asarray([ang(patch_normal(t), g) for t, g in zip(ts, n_gt)])
        pe = np.linalg.norm(X - X_gt, axis=1) / half
        z = [median_zncc(t) for t in ts]
        line = (
            f"{label:9s} normal err med {np.median(ne):5.1f} p90 "
            f"{np.percentile(ne, 90):5.1f} | pos err (halves) med {np.median(pe):.3f} "
            f"p90 {np.percentile(pe, 90):.3f} | zncc med {np.median(z):.3f}"
        )
        if sloped.any():
            line += f" | sloped normal err med {np.median(ne[sloped]):5.1f}"
        print(line)

    def neighbour_normal(i, X, current):
        r = args.nb_halves * half[i]
        near = [j for j in range(len(X)) if np.linalg.norm(X[j] - X[i]) < r]
        if len(near) < 3:
            return current
        P = X[near] - X[near].mean(0)
        vals, vecs = np.linalg.eigh(P.T @ P / len(near))
        # The in-plane directions the positions pin down.
        fixed = [
            vecs[:, k] for k in range(3) if np.sqrt(max(vals[k], 0)) > 0.25 * half[i]
        ]
        if len(fixed) >= 2:
            n = unit(np.cross(fixed[-1], fixed[-2]))
        elif len(fixed) == 1:
            a = fixed[0]
            n = unit(current - (current @ a) * a)
        else:
            return current
        return n if n @ current > 0 else -n

    ts = [tilt(t, init_normal(i, t)) for i, t in enumerate(tracks)]
    report("start", ts)
    ts = [fit(t) for t in ts]
    report("fit", ts)
    for rnd in range(1, args.rounds + 1):
        X = np.asarray([t.position for t in ts])
        ns = [neighbour_normal(i, X, unit(patch_normal(t))) for i, t in enumerate(ts)]
        ts = [tilt(t, n) for t, n in zip(ts, ns)]
        report(f"tilt {rnd}", ts)
        ts = [fit(t) for t in ts]
        report(f"refit {rnd}", ts)


if __name__ == "__main__":
    main()
