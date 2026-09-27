# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""A group of ground-truth points around one: views, readings, photometry, and run outcomes.

For each finite point within ``--radius`` metres of ``CENTRE``: its view count
and widest baseline, the per-view ZNCC of the ground-truth track, its normal's
angle to the group's mean normal, the photometric normal's error (with
``renormal``'s settings, at the ground truth's position and views), and the
pose on a coarse grid of tilts and depths with the best median ZNCC, against
the ground truth's own. With ``--runs``, each run's good tracks over queries
and median normal error on the point, in each pass.

    pixi run -e test python scripts/track_at_pixel/surface/group.py DATASET CENTRE [--radius 1.5] [--runs RUN ...]
"""

import argparse
import json
from pathlib import Path

import numpy as np

from _common import ang, load, tangent, unit, zncc_and_shift


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dataset")
    ap.add_argument("centre", type=int)
    ap.add_argument("--radius", type=float, default=1.5)
    ap.add_argument("--runs", nargs="*", default=[], type=Path)
    args = ap.parse_args()

    from sfmtool._sfmtool import bench as B

    from candidates.common import patch_normal
    from candidates.renormal import camera_views

    ds, edited = load(args.dataset)
    c = ds.point_xyz[args.centre]
    group = [
        int(p)
        for p in np.flatnonzero(ds.point_w != 0)
        if np.linalg.norm(ds.point_xyz[p] - c) < args.radius
    ]
    group.sort(key=lambda p: np.linalg.norm(ds.point_xyz[p] - c))

    runs = {}
    for run in args.runs:
        rows = {}
        for line in open(run / "rows.jsonl"):
            r = json.loads(line)
            if r["point"] in group:
                rows.setdefault((r["point"], r["pass"]), []).append(r)
        runs[run.name] = rows

    def pose(gt, cam, X0, normal, factor):
        t, _ = B.tilt_patch(gt, edited, tuple(float(x) for x in normal))
        pl = t.placement
        u, v = unit(pl["u_halfvec"]), unit(pl["v_halfvec"])
        n = unit(np.cross(u, v))
        d = (factor - 1.0) * (X0 - cam.center)
        t, _ = B.translate_patch(t, edited, (float(d @ u), float(d @ v), float(d @ n)))
        return t

    def photo(X, imgs, half, init):
        from sfmtool._sfmtool import patches

        e1, e2 = tangent(init)
        cloud = patches.PatchCloud.from_halfvec_arrays(
            np.asarray([e1 * half], np.float32),
            np.asarray([e2 * half], np.float32),
            np.asarray([X], np.float64),
        )
        r = cloud.refine_normals(
            camera_views(ds.holdout(None)),
            ds.pyramids,
            view_indices=[imgs],
            use_stored_keypoints=False,
            min_views=2,
            angular_range_deg=45.0,
            fronto_prior_weight=0.05,
        )
        return unit(r["normal"][0])

    ref = ds.point_normal[args.centre]
    mean_n = unit(
        sum(ds.point_normal[p] * np.sign(ds.point_normal[p] @ ref) for p in group)
    )
    print(f"{len(group)} points within {args.radius} m of {args.centre}\n")
    head = "pt   dist  views maxbase  GT zncc per view            vs-grp photo  best grid (f, tilt@az) z   GT z"
    for name in runs:
        head += f" | {name}: full, empty good/queries med normal err"
    print(head)
    for p in group:
        _, gt = B.create_track(B.Bench(), edited, p)
        imgs = [int(o["image"]) for o in gt.observations]
        X0 = ds.point_xyz[p]
        n_gt = unit(patch_normal(gt))
        cam = ds.cameras[imgs[0]]
        if n_gt @ (cam.center - X0) < 0:
            n_gt = -n_gt
        dirs = [unit(ds.cameras[i].center - X0) for i in imgs]
        maxbase = max(ang(a, b) for a in dirs for b in dirs)
        z, _ = zncc_and_shift(ds, edited, gt)
        try:
            ph_err = ang(
                photo(X0, imgs, 2 * float(ds.point_half[p]), unit(sum(dirs))), n_gt
            )
        except (ValueError, RuntimeError):
            ph_err = float("nan")
        e1, e2 = tangent(n_gt)
        best = (-1.0, (1.0, 0, 0))
        for tilt in (0, 10, 20, 30):
            for az in [0] if tilt == 0 else range(0, 360, 60):
                r_, th = np.radians(tilt), np.radians(az)
                n = unit(
                    np.cos(r_) * n_gt + np.sin(r_) * (np.cos(th) * e1 + np.sin(th) * e2)
                )
                for f in (0.97, 0.985, 1.0, 1.015, 1.03):
                    try:
                        zz, _ = zncc_and_shift(ds, edited, pose(gt, cam, X0, n, f))
                    except ValueError:
                        continue
                    m = float(np.nanmedian(zz))
                    if m > best[0]:
                        best = (m, (f, tilt, az))
        cells = []
        for name in runs:
            parts = []
            for pas in ("full", "empty"):
                rows = runs[name].get((p, pas), [])
                good = [r for r in rows if r.get("good")]
                ne = [
                    r["normal_err_deg"]
                    for r in good
                    if r.get("normal_err_deg") is not None
                ]
                med = np.median(ne) if ne else float("nan")
                parts.append(f"{len(good)}/{len(rows)} {med:4.0f}")
            cells.append("  ".join(parts))
        f, tilt, az = best[1]
        print(
            f"{p:4d} {np.linalg.norm(X0 - c):5.2f} {len(imgs):4d} {maxbase:7.1f}  "
            f"{' '.join(f'{x:.2f}' for x in z):28s} {ang(n_gt, mean_n):5.1f} {ph_err:5.1f}  "
            f"{f:.3f},{tilt:2d}@{az:3d} {best[0]:.2f}  {np.nanmedian(z):.2f}"
            + "".join(f" | {c}" for c in cells)
        )


if __name__ == "__main__":
    main()
