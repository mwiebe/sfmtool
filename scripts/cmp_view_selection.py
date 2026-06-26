#!/usr/bin/env python3
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Head-to-head: prototype ``vet()`` (Python) vs production ``select_views`` (Rust).

Validates that the shipped view-selection kernel reproduces the prototype's
photometric vetting. On the same reconstruction, same points, same (track-)
refined normals and patch frames, both paths score every geometrically-visible
candidate view against a track-seeded reference and admit those clearing the
adaptive threshold. Reports, per point, how many extra views each admits and the
Jaccard overlap of the two admitted sets (1.00 = identical decisions).

Usage::

    pixi run python scripts/cmp_view_selection.py RECON.sfmr [-n POINTS]
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
    vet,
)

from sfmtool._solve_strips import _SolveStrips  # noqa: E402
from sfmtool._sfmtool import SfmrReconstruction  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("recon", type=Path, help="solved .sfmr")
    ap.add_argument("-n", "--num", type=int, default=12, help="points to compare")
    args = ap.parse_args()

    recon = SfmrReconstruction.load(str(args.recon))
    s = _SolveStrips(recon, recon.workspace_dir, patch=PATCH, extent_factor=EXTENT_FACTOR)
    imgs = [s.image(i) for i in range(len(s.names))]
    w = gauss_window(PATCH)
    cidx = {int(p): k for k, p in enumerate(s.cloud.point_ids)}
    finite = np.abs(np.asarray(recon.positions_xyzw)[:, 3]) > 1e-9

    # Candidate points: finite, well-tracked, with extra geometrically-visible
    # views to vet (the points where selection has a decision to make).
    cand = []
    for pid in s.cloud.point_ids:
        pid = int(pid)
        if not finite[pid]:
            continue
        tracked = sorted(set(s.obs.get(pid, [])))
        if len(tracked) < 3:
            continue
        vis = geometric_views(s, pid)
        extra = len(set(vis) - set(tracked))
        if extra:
            cand.append((pid, set(tracked), vis, extra))
    if not cand:
        raise SystemExit("no finite, well-tracked points with extra visible views")
    cand.sort(key=lambda t: t[3], reverse=True)
    pool = cand[: args.num]
    ids = [pid for pid, *_ in pool]
    print(f"{args.recon}")
    print(f"points: {len(ids)} (finite, track>=3, with extra visible views)\n")

    # Shared: refine the SAME cloud's normals over the track so both paths vet
    # under identical normals/frames. Production resolution set to PATCH to match
    # the prototype's PATCH-sized scored core.
    s.cloud.refine_normals(
        recon, imgs, point_ids=ids, resolution=PATCH,
        angular_range_deg=ANG_RANGE, init_steps=INIT_STEPS,
    )

    # Production: select_views (Rust), same threshold and resolution.
    sels = {
        int(d["point_id"]): d
        for d in s.cloud.select_views(
            recon, imgs, min_relative_zncc=0.7, resolution=PATCH, point_ids=ids
        )
    }

    hdr = (f"{'pid':>5} {'trk':>3} {'cand':>4} | {'proto+':>6} {'prod+':>5} "
           f"{'J(extra)':>8} | {'base':>5} {'self':>5}")
    print(hdr)
    print("-" * len(hdr))
    jacc, dproto, dprod = [], [], []
    for pid, T, vis, extra in pool:
        normal = np.asarray(s.cloud[cidx[pid]].normal, float)
        # Prototype vet: admit extra views whose corr >= 0.7 * track self-corr.
        _, _, _, admitted_p, base = vet(
            s, pid, normal, vis, T, sorted(T), w, PATCH, s._half(pid), 0.7, 0.1
        )
        proto_extra = {int(i) for i in admitted_p}  # vet returns extras only
        d = sels.get(pid)
        prod_all = {int(i) for i in np.asarray(d["admitted"]).tolist()}
        prod_extra = prod_all - set(T)
        self_agree = float(d["self_agreement"])
        union = len(proto_extra | prod_extra)
        j = len(proto_extra & prod_extra) / union if union else 1.0
        jacc.append(j)
        dproto.append(len(proto_extra))
        dprod.append(len(prod_extra))
        print(
            f"{pid:>5} {len(T):>3} {extra:>4} | {len(proto_extra):>6} "
            f"{len(prod_extra):>5} {j:>8.2f} | {base:>5.2f} {self_agree:>5.2f}"
        )

    print("-" * len(hdr))
    print(f"mean admitted-extra:  prototype {np.mean(dproto):.1f}   "
          f"production {np.mean(dprod):.1f}")
    print(f"mean Jaccard(extra admitted sets): {np.mean(jacc):.2f}")


if __name__ == "__main__":
    main()
