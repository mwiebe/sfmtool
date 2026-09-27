# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""The depth x normal landscape around one ground-truth point.

The ground-truth track is tilted from its own normal by 5 to 30 degrees in
four directions, and moved along the queried view's ray to 0.97 to 1.03 times
its depth. At each pose the sightings are the patch's projections (the view
set is held fixed), and the table gives the median ZNCC and the median
correlation peak offset in pixels. The last row starts from the mean viewing
direction.

    pixi run -e test python scripts/track_at_pixel/surface/landscape.py DATASET POINT [QUERY_OBS]
"""

import argparse

import numpy as np

from _common import ang, load, tangent, unit, zncc_and_shift


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dataset")
    ap.add_argument("point", type=int)
    ap.add_argument("query_obs", type=int, nargs="?", default=0)
    args = ap.parse_args()

    from sfmtool._sfmtool import bench as B

    from candidates.common import patch_normal

    ds, edited = load(args.dataset)
    _, gt = B.create_track(B.Bench(), edited, args.point)
    img_q = int(gt.observations[args.query_obs]["image"])
    cam = ds.cameras[img_q]
    X0 = np.asarray(gt.position, float)
    n_gt = unit(patch_normal(gt))
    if n_gt @ (cam.center - X0) < 0:
        n_gt = -n_gt
    views = [int(o["image"]) for o in gt.observations]
    mv = unit(sum(unit(ds.cameras[i].center - X0) for i in views))
    print(f"point {args.point}: {len(views)} views {views}, query image {img_q}")
    print(f"GT normal {np.round(n_gt, 3)}; mean view {ang(n_gt, mv):.1f} deg off GT")
    print(
        f"half size {np.linalg.norm(gt.placement['u_halfvec']):.3f} m, "
        f"depth {cam.depth(X0):.2f} m"
    )
    for i in views:
        c = ds.cameras[i]
        print(
            f"  view {i}: angle to GT normal {ang(n_gt, c.center - X0):5.1f}, "
            f"dist {np.linalg.norm(c.center - X0):.2f}"
        )

    def pose(normal, factor):
        t, _ = B.tilt_patch(gt, edited, tuple(float(x) for x in normal))
        pl = t.placement
        u, v = unit(pl["u_halfvec"]), unit(pl["v_halfvec"])
        n = unit(np.cross(u, v))
        d = (factor - 1.0) * (X0 - cam.center)
        t, _ = B.translate_patch(t, edited, (float(d @ u), float(d @ v), float(d @ n)))
        return t

    z, s = zncc_and_shift(ds, edited, gt)
    print("\nGT as stored: zncc", np.round(z, 3), "peak shift px", np.round(s, 2))

    e1, e2 = tangent(n_gt)
    factors = [0.97, 0.98, 0.99, 0.995, 1.0, 1.005, 1.01, 1.02, 1.03]

    def row(n):
        cells = []
        for f in factors:
            try:
                z, s = zncc_and_shift(ds, edited, pose(n, f))
                cells.append(f"{np.nanmedian(z):5.2f}/{np.nanmedian(s):4.1f}")
            except ValueError:
                cells.append("     --    ")
        return " ".join(cells)

    print("\nmedian ZNCC / median peak shift (px) at the projections, view set fixed")
    print("tilt    " + " ".join(f"{f:>10.3f}" for f in factors))
    for tilt in (0, 5, 10, 20, 30):
        for az in [0] if tilt == 0 else range(0, 360, 90):
            r, th = np.radians(tilt), np.radians(az)
            n = unit(np.cos(r) * n_gt + np.sin(r) * (np.cos(th) * e1 + np.sin(th) * e2))
            print(f"{tilt:2d}@{az:3d}  " + row(n))
    print("meanview " + row(mv))


if __name__ == "__main__":
    main()
