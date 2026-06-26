#!/usr/bin/env python3
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Head-to-head: prototype ``congeal()`` (Python) vs production ``localize_keypoints``
(Rust).

Validates the shipped keypoint-localization kernel against the prototype's
congealing loop (group-wise sub-pixel translation registration). On the same
points, same (track-)refined normals, and the same admitted view set fed to both,
reports the median total sub-pixel shift each path finds, the before/after
leave-one-out ZNCC, and how many views each keeps.

The prototype probe congeals *all* given views; production additionally drops
views that won't co-register in-loop (grazing, ``max_shift_px``, low LOO), so it
typically reaches a slightly higher final LOO over fewer, cleaner views.

Usage::

    pixi run python scripts/cmp_keypoint_localization.py RECON.sfmr [-n POINTS]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _normal_strip_lib import (  # noqa: E402
    ANG_RANGE,
    EXTENT_FACTOR,
    INIT_STEPS,
    PATCH,
    gauss_window,
    geometric_views,
)
from exp_reference_refine import congeal  # noqa: E402

from sfmtool._solve_strips import _SolveStrips  # noqa: E402
from sfmtool._sfmtool import SfmrReconstruction  # noqa: E402

RENDER_RES = 64
OFF = (RENDER_RES - PATCH) // 2
SEARCH = 6
ITERS = 5


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("recon", type=Path, help="solved .sfmr")
    ap.add_argument("-n", "--num", type=int, default=10, help="points to compare")
    args = ap.parse_args()

    recon = SfmrReconstruction.load(str(args.recon))
    s = _SolveStrips(recon, recon.workspace_dir, patch=PATCH, extent_factor=EXTENT_FACTOR)
    imgs = [s.image(i) for i in range(len(s.names))]
    w = gauss_window(PATCH)
    cidx = {int(p): k for k, p in enumerate(s.cloud.point_ids)}
    finite = np.abs(np.asarray(recon.positions_xyzw)[:, 3]) > 1e-9

    cand = []
    for pid in s.cloud.point_ids:
        pid = int(pid)
        if not finite[pid]:
            continue
        tracked = sorted(set(s.obs.get(pid, [])))
        if len(tracked) < 3:
            continue
        vis = geometric_views(s, pid)
        if len(set(vis) - set(tracked)):
            cand.append((pid, set(tracked), vis, len(set(vis) - set(tracked))))
    if not cand:
        raise SystemExit("no finite, well-tracked points with extra visible views")
    cand.sort(key=lambda t: t[3], reverse=True)
    pool = cand[: args.num]
    ids = [pid for pid, *_ in pool]

    # Shared: refine normals, then pick a common admitted view set per point via
    # production select_views, so both congealers operate on the same stack.
    s.cloud.refine_normals(
        recon, imgs, point_ids=ids, resolution=PATCH,
        angular_range_deg=ANG_RANGE, init_steps=INIT_STEPS,
    )
    sels = {
        int(d["point_id"]): np.asarray(d["admitted"]).tolist()
        for d in s.cloud.select_views(
            recon, imgs, min_relative_zncc=0.7, resolution=PATCH, point_ids=ids
        )
    }

    # Production congealing over those exact view sets.
    loc = {
        int(d["point_id"]): d
        for d in s.cloud.localize_keypoints(
            recon, imgs, view_sets=sels, max_iters=ITERS, search=float(SEARCH),
            max_shift_px=3.0, min_relative_zncc=0.7, resolution=PATCH,
        )
    }

    print(f"{args.recon}")
    print(f"points: {len(ids)}   (proto keeps all views; prod may drop)\n")
    hdr = (f"{'pid':>5} {'G':>3} | {'proto shift':>11} {'prod shift':>10} | "
           f"{'proto LOO 0>N':>14} {'prod LOO N':>10} | {'prod kept':>9}")
    print(hdr)
    print("-" * len(hdr))
    ps, qs, pl, ql = [], [], [], []
    for pid, T, vis, extra in pool:
        G = sorted(int(i) for i in sels[pid])
        if len(G) < 2:
            continue
        normal = np.asarray(s.cloud[cidx[pid]].normal, float)
        ext = s._half(pid) * (RENDER_RES / PATCH)
        stats, _, _ = congeal(
            s, pid, normal, G, ext, RENDER_RES, OFF, w,
            iters=ITERS, search=SEARCH, coarse=True,
        )
        proto_shift = stats["shift_src"]
        proto_loo0, proto_looN = stats["loo"]
        d = loc.get(pid)
        off_px = np.asarray(d["offsets_px"], float)
        prod_shift = float(np.median(off_px)) if len(off_px) else float("nan")
        prod_loo = float(np.nanmean(np.asarray(d["loo_zncc"], float)))
        kept = len(d["views"])
        ps.append(proto_shift)
        qs.append(prod_shift)
        pl.append(proto_looN)
        ql.append(prod_loo)
        print(
            f"{pid:>5} {len(G):>3} | {proto_shift:>10.2f}p {prod_shift:>9.2f}p | "
            f"{proto_loo0:>6.3f}>{proto_looN:>6.3f} {prod_loo:>10.3f} | "
            f"{kept:>4}/{len(G):<4}"
        )

    print("-" * len(hdr))
    print(f"median shift (src px):   prototype {np.median(ps):.2f}   "
          f"production {np.median(qs):.2f}")
    print(f"mean final LOO ZNCC:     prototype {np.mean(pl):.3f}   "
          f"production {np.mean(ql):.3f}")


if __name__ == "__main__":
    main()
