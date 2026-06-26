#!/usr/bin/env python3
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment: MVS-style normal refinement over all geometrically-visible views.

Normal refinement (``sfm xform --refine-normals``) scores a point's oriented
surfel over the images that *matched a feature* there — the point's track. This
experiment instead feeds the refiner **every image that geometrically sees the
point** (projects in front of the camera and inside the frame), turning the
sparse-feature consensus into a denser, MVS-like one. We then render the result
as padded patch strips so we can see two things at a glance:

* **how the normal refines more broadly** — the consensus normal found over the
  full visible set vs. the track-only set (the header reports the angle between
  them and the consensus photoconsistency each achieves over the visible set); and
* **how the robust inclusion weights work out** — each view tile is annotated
  with its IRLS (Tukey) consensus weight, so the views the robust estimator
  trusts vs. down-weights (occluded / wrong-surface / off-plane) are visible.

See ``exp_goodview_normal_strips.py`` for the staged alternative (refine on the
track first, then grow the view set by agreement) that this experiment motivates.

Usage::

    pixi run python scripts/exp_mvs_normal_strips.py RECON.sfmr -o out.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import cv2
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _normal_strip_lib import (  # noqa: E402
    ANG_RANGE,
    EXTENT_FACTOR,
    INIT_STEPS,
    PATCH,
    gauss_window,
    geometric_views,
    irls_consensus,
    legend_bar,
    render_cores,
    strip_image,
)

from sfmtool._solve_strips import _SolveStrips  # noqa: E402
from sfmtool._sfmtool import PatchCloud, SfmrReconstruction  # noqa: E402

TRACK_COL = (0, 255, 255)  # yellow: a feature-matched (track) view
EXTRA_COL = (80, 80, 255)  # red: visible but not matched (MVS-only)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("recon", type=Path)
    ap.add_argument("-o", "--output", type=Path, default=Path("mvs_normal_strips.png"))
    ap.add_argument("-n", "--num", type=int, default=8, help="points to render")
    ap.add_argument("--context", type=int, default=96, help="padded render resolution (px)")
    args = ap.parse_args()

    recon = SfmrReconstruction.load(str(args.recon))
    s = _SolveStrips(recon, recon.workspace_dir, patch=PATCH, extent_factor=EXTENT_FACTOR)
    imgs = [s.image(i) for i in range(len(s.names))]
    render_res = max(args.context, PATCH)
    w = gauss_window(PATCH)

    finite = np.abs(np.asarray(recon.positions_xyzw)[:, 3]) > 1e-9

    # Candidate points: finite, well-tracked, and seen by extra (non-track) views,
    # ranked by how much the MVS expansion adds (most extra coverage first).
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
        if extra > 0:
            cand.append((pid, tracked, vis, extra))
    cand.sort(key=lambda t: t[3], reverse=True)
    sel = cand[: args.num]
    if not sel:
        raise SystemExit("no finite, well-tracked points with extra visible views")

    print(f"rendering {len(sel)} points (of {len(cand)} with extra coverage)")
    sel_ids = [pid for pid, *_ in sel]
    vis_of = {pid: vis for pid, _, vis, _ in sel}
    track_of = {pid: set(tr) for pid, tr, _, _ in sel}

    # Per-patch view-index override: the geometric (MVS) view set for selected
    # points, empty elsewhere (point_ids restricts which patches actually refine).
    cidx = {int(p): k for k, p in enumerate(s.cloud.point_ids)}
    view_idx = [[] for _ in s.cloud.point_ids]
    for pid in sel_ids:
        view_idx[cidx[pid]] = [int(i) for i in vis_of[pid]]

    common = dict(resolution=PATCH, angular_range_deg=ANG_RANGE, init_steps=INIT_STEPS)

    # Track-only refinement (default view set).
    cloud_t = PatchCloud.from_reconstruction(recon, normal="stored", extent_value=EXTENT_FACTOR)
    res_t = cloud_t.refine_normals(recon, imgs, point_ids=sel_ids, **common)
    ci_t = {int(p): k for k, p in enumerate(cloud_t.point_ids)}
    n_tracks = {pid: np.asarray(res_t["normal"])[ci_t[pid]] for pid in sel_ids}

    # MVS refinement over all geometrically-visible views.
    cloud_m = PatchCloud.from_reconstruction(recon, normal="stored", extent_value=EXTENT_FACTOR)
    res_m = cloud_m.refine_normals(
        recon, imgs, point_ids=sel_ids, view_indices=view_idx, **common
    )
    n_mvs = {pid: np.asarray(res_m["normal"])[cidx[pid]] for pid in sel_ids}

    blocks: list[np.ndarray] = []
    for pid in sel_ids:
        vis = vis_of[pid]
        ext = s._half(pid) * (render_res / PATCH)
        # Render the same MVS view set under both normals so the consensus
        # tightening from track-only -> MVS is directly comparable.
        tiles_t, cores_t = render_cores(s, pid, n_tracks[pid], vis, render_res, ext)
        tiles_m, cores_m = render_cores(s, pid, n_mvs[pid], vis, render_res, ext)
        wt_t, phi_t = irls_consensus(cores_t, w)
        wt_m, phi_m = irls_consensus(cores_m, w)
        dang = float(np.degrees(np.arccos(np.clip(
            np.dot(n_tracks[pid] / np.linalg.norm(n_tracks[pid]),
                   n_mvs[pid] / np.linalg.norm(n_mvs[pid])), -1, 1))))

        cols = [TRACK_COL if i in track_of[pid] else EXTRA_COL for i in vis]
        v = len(vis)
        rel_t = [float(np.clip(x * v, 0, 2) / 2) for x in wt_t]
        rel_m = [float(np.clip(x * v, 0, 2) / 2) for x in wt_m]
        row_t = strip_image(tiles_t, vis, cols, [f"{x:.2f}" for x in wt_t], rel_t,
                            render_res, label="track normal")
        row_m = strip_image(tiles_m, vis, cols, [f"{x:.2f}" for x in wt_m], rel_m,
                            render_res, label="MVS normal")

        header = np.full((30, row_t.shape[1], 3), 18, np.uint8)
        msg = (f"pt {pid}  views: {len(track_of[pid])} track + "
               f"{v - len(track_of[pid])} extra = {v} MVS   "
               f"normal moved {dang:.1f} deg   "
               f"consensus Phi: track-n={phi_t:+.3f} -> MVS-n={phi_m:+.3f}")
        cv2.putText(header, msg, (6, 20), cv2.FONT_HERSHEY_SIMPLEX, 0.46, (180, 220, 180), 1, cv2.LINE_AA)
        blocks.extend([header, row_t, row_m, np.full((6, header.shape[1], 3), 18, np.uint8)])

        uni = 1.0 / v
        print(f"pt {pid:4d}: {len(track_of[pid]):2d} track +{v-len(track_of[pid]):2d} "
              f"extra | normal {dang:5.1f}deg | Phi {phi_t:+.3f}->{phi_m:+.3f} | "
              f"w[min={wt_m.min():.3f} max={wt_m.max():.3f} uni={uni:.3f}]")

    width = max(b.shape[1] for b in blocks)
    blocks = [np.pad(b, ((0, 0), (0, width - b.shape[1]), (0, 0)), constant_values=18) for b in blocks]
    legend = legend_bar(width, [
        (TRACK_COL, "track view (matched)"),
        (EXTRA_COL, "extra view (visible, not matched)"),
        ((255, 255, 255), "validated extent"),
    ])
    montage = np.vstack([legend, *blocks])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(args.output), montage)
    print(f"wrote {args.output} ({montage.shape[1]}x{montage.shape[0]})")


if __name__ == "__main__":
    main()
