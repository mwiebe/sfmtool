# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Two half-size patches split along an in-plane axis read the tilt about the other.

For each ground-truth point, the normal starts at the ground truth's tilted by
-15, 0 and +15 degrees about the patch's horizontal in-plane axis (so it leans
up or down). Each round the patch is tilted to the current normal and fitted,
then two copies of half the size are moved half a size up and down in its
plane (along the in-plane direction closest to world up) and fitted, and the
line between their fitted centres sets the tilt: the normal is made
perpendicular to it, keeping its horizontal part. Reported per round: the
vertical-lean error against the ground truth, and the total normal error.

    pixi run -e test python scripts/track_at_pixel/surface/split_patch.py DATASET POINT [POINT ...]
"""

import argparse

import numpy as np

from _common import UP, ang, load, move, tilt, unit


def lean(n):
    """Signed angle of the normal above the horizontal plane."""
    return float(np.degrees(np.arcsin(np.clip(unit(n) @ UP, -1, 1))))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dataset")
    ap.add_argument("points", type=int, nargs="+")
    args = ap.parse_args()

    from sfmtool._sfmtool import bench as B

    from candidates.common import patch_normal

    ds, edited = load(args.dataset)

    def split_normal(t, n, half):
        up_in_plane = unit(UP - (UP @ n) * n)
        centres = []
        for sign in (1, -1):
            s, _ = B.resize_patch(t, edited, half / 2)
            s = move(edited, s, sign * (half / 2) * up_in_plane)
            s, _ = B.fit(s, edited, ds.pyramids)
            centres.append(np.asarray(s.position, float))
        a = unit(centres[0] - centres[1])
        out = unit(n - (n @ a) * a)
        return out if out @ n > 0 else -out

    for p in args.points:
        _, gt = B.create_track(B.Bench(), edited, p)
        n_gt = unit(patch_normal(gt))
        cam = ds.cameras[int(gt.observations[0]["image"])]
        if n_gt @ (cam.center - np.asarray(gt.position)) < 0:
            n_gt = -n_gt
        half = float(np.linalg.norm(gt.placement["u_halfvec"]))
        horiz = unit(np.cross(UP, n_gt))
        print(f"point {p}: half {half:.3f} m, GT lean {lean(n_gt):+.1f} deg")
        for start in (-15.0, 0.0, 15.0):
            r = np.radians(start)
            n = unit(np.cos(r) * n_gt + np.sin(r) * np.cross(horiz, n_gt))
            t = gt
            trace = [f"{lean(n) - lean(n_gt):+5.1f}"]
            try:
                for _ in range(3):
                    t = tilt(ds, edited, t, n)
                    t, _ = B.fit(t, edited, ds.pyramids)
                    n = split_normal(t, n, half)
                    trace.append(f"{lean(n) - lean(n_gt):+5.1f}")
            except ValueError as e:
                trace.append(f"refused: {e}"[:60])
            print(
                f"  start {start:+5.1f}: lean error by round {' -> '.join(trace)}; "
                f"total error {ang(n, n_gt):.1f}"
            )


if __name__ == "__main__":
    main()
