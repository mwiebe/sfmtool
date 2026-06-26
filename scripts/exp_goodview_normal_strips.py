#!/usr/bin/env python3
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment: staged good-view normal refinement.

The companion ``exp_mvs_normal_strips.py`` showed that throwing *every*
geometrically-visible view at the refiner admits self-occluded views the robust
consensus can't reject (its MAD scale inflates and the Tukey cutoff stops
discriminating). This experiment does the staged thing instead — the spec's
iterative good-view set (``specs/core/patch-normal-refinement.md`` item 2):

1. **Refine on the track.** The matched views give a trustworthy seed normal and
   a trustworthy appearance template.
2. **Score the rest against it.** Render every other geometrically-visible view
   under that normal and correlate it (windowed ZNCC) to the track template.
3. **Admit the good candidates.** Accept the extra views whose correlation clears
   a threshold tied to the track's own self-agreement (so the bar adapts to how
   textured/consistent the surfel is) — these are the genuinely-visible,
   non-occluded views.
4. **Refine again** over track + accepted, and repeat a couple of rounds (each
   round can admit more views as the template sharpens).

Rendered as padded patch strips: the top strip is the track normal (each extra
view annotated with its correlation to the track template, green = admitted, red
= rejected/occluded); the bottom strip is the final good-view normal (each
admitted view annotated with its robust inclusion weight).

Usage::

    pixi run python scripts/exp_goodview_normal_strips.py RECON.sfmr -o out.png
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
    _phi,
    gauss_window,
    geometric_views,
    irls_weights,
    legend_bar,
    render_cores,
    strip_image,
    vet,
    znorm_stack,
)

from sfmtool._solve_strips import _SolveStrips  # noqa: E402
from sfmtool._sfmtool import PatchCloud, SfmrReconstruction  # noqa: E402

TRACK_COL = (0, 255, 255)  # yellow: feature-matched (track) view
GOOD_COL = (90, 230, 90)  # green: admitted by photometric agreement
DROP_COL = (80, 80, 255)  # red: rejected (occluded / wrong surface)


def _angle(a, b) -> float:
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    return float(np.degrees(np.arccos(np.clip(np.dot(a, b), -1.0, 1.0))))


def _report_sift_overlap(s, recon, sel_ids, track_of, good) -> None:
    """How many admitted (added) views actually have a detected SIFT keypoint near
    the point's projection? Distinguishes views the matcher merely *missed* (a
    feature is there) from genuinely feature-less coverage the photometric vetting
    uniquely recovers. Track views (themselves keypoints) are the ~0px baseline."""
    from sfmtool.sift.file import SiftReader, get_sift_path_from_recon

    kp_cache: dict[int, np.ndarray] = {}

    def kps(i):
        if i not in kp_cache:
            path = get_sift_path_from_recon(recon, s.names[i])
            kp_cache[i] = np.asarray(SiftReader(str(path)).read_positions(), float)
        return kp_cache[i]

    def nearest(i, center):
        pc = s.rot_of[i] @ (center - s.centers[i])
        if pc[2] <= 1e-6:
            return None
        u, v = s.cam_of[i].project(pc[0] / pc[2], pc[1] / pc[2])
        k = kps(i)
        if len(k) == 0:
            return None
        return float(np.min(np.hypot(k[:, 0] - u, k[:, 1] - v)))

    adm, trk = [], []
    for pid in sel_ids:
        center = s.positions[int(pid)]
        for i in sorted(set(good[pid]) - track_of[pid]):
            d = nearest(i, center)
            if d is not None:
                adm.append(d)
        for i in sorted(track_of[pid]):
            d = nearest(i, center)
            if d is not None:
                trk.append(d)
    adm, trk = np.asarray(adm), np.asarray(trk)
    if not len(adm):
        print("\nSIFT overlap: no admitted views")
        return
    print(f"\nSIFT overlap of {len(adm)} ADMITTED (added) views "
          f"[median nearest-keypoint {np.median(adm):.1f}px]:")
    for t in (2.0, 5.0, 10.0):
        print(f"    within {t:>2.0f}px of a detected keypoint: "
              f"{int((adm <= t).sum()):3d}/{len(adm)} ({100 * (adm <= t).mean():.0f}%)")
    if len(trk):
        print(f"  baseline — track views (are keypoints): median {np.median(trk):.1f}px, "
              f"{100 * (trk <= 2).mean():.0f}% within 2px")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("recon", type=Path)
    ap.add_argument("-o", "--output", type=Path, default=Path("goodview_normal_strips.png"))
    ap.add_argument("-n", "--num", type=int, default=8, help="points to render")
    ap.add_argument("--context", type=int, default=96, help="padded render resolution (px)")
    ap.add_argument("--iters", type=int, default=2, help="grow/refine rounds")
    ap.add_argument("--accept-frac", type=float, default=0.7,
                    help="admit a view if its corr >= accept_frac * mean track self-corr")
    ap.add_argument("--floor", type=float, default=0.1, help="absolute corr floor to admit")
    ap.add_argument("--rank", choices=("admits", "extra"), default="admits",
                    help="select points to render by admitted-view count (default) or "
                    "by raw extra-view coverage")
    ap.add_argument("--pool", type=int, default=200,
                    help="candidate points pre-vetted when --rank admits")
    ap.add_argument("--max-views", type=int, default=0,
                    help="cap tiles shown per strip (evenly subsampled; 0 = show all). "
                    "Decisions/Phi are always over the full visible set; this only "
                    "keeps wide many-view montages legible.")
    ap.add_argument("--sift-overlap", action="store_true",
                    help="report how many admitted (added) views have a detected SIFT "
                    "keypoint near the point's projection — i.e. views the matcher "
                    "missed vs. genuinely feature-less coverage.")
    args = ap.parse_args()

    recon = SfmrReconstruction.load(str(args.recon))
    s = _SolveStrips(recon, recon.workspace_dir, patch=PATCH, extent_factor=EXTENT_FACTOR)
    imgs = [s.image(i) for i in range(len(s.names))]
    render_res = max(args.context, PATCH)
    w = gauss_window(PATCH)
    common = dict(resolution=PATCH, angular_range_deg=ANG_RANGE, init_steps=INIT_STEPS)

    finite = np.abs(np.asarray(recon.positions_xyzw)[:, 3]) > 1e-9
    cidx = {int(p): k for k, p in enumerate(s.cloud.point_ids)}

    # Candidate pool: finite, well-tracked, with extra geometrically-visible views
    # to vet. Capped (by extra coverage) to bound the pre-vet cost.
    cand = []
    for pid in s.cloud.point_ids:
        pid = int(pid)
        if not finite[pid]:
            continue
        tracked = sorted(set(s.obs.get(pid, [])))
        if len(tracked) < 3:
            continue
        vis = geometric_views(s, pid)
        if set(vis) - set(tracked):
            cand.append((pid, set(tracked), vis, len(set(vis) - set(tracked))))
    if not cand:
        raise SystemExit("no finite, well-tracked points with extra visible views")
    cand.sort(key=lambda t: t[3], reverse=True)
    pool = cand[: args.pool] if args.rank == "admits" else cand[: args.num]
    pool_ids = [pid for pid, *_ in pool]
    ext_of = {pid: s._half(pid) * (render_res / PATCH) for pid in pool_ids}

    # Stage 1: track-only refinement (default view set) — the trusted seed, for the
    # whole pool so we can rank by how many extra views each point ends up admitting.
    cloud_t = PatchCloud.from_reconstruction(recon, normal="stored", extent_value=EXTENT_FACTOR)
    res_t = cloud_t.refine_normals(recon, imgs, point_ids=pool_ids, **common)
    n_track = {pid: np.asarray(res_t["normal"])[cidx[pid]] for pid in pool_ids}

    if args.rank == "admits":
        # Pre-vet the pool (round 0) and keep the points that admit the most extra
        # views — the cases that best show the expansion working, not just the
        # heaviest-occlusion cases that raw coverage surfaces.
        scored = []
        for pid, T, vis, extra in pool:
            _, _, _, admitted, _ = vet(
                s, pid, n_track[pid], vis, T, sorted(T), w,
                render_res, ext_of[pid], args.accept_frac, args.floor,
            )
            scored.append((pid, T, vis, extra, len(admitted)))
        scored.sort(key=lambda r: (r[4], r[3]), reverse=True)
        sel = scored[: args.num]
        print(f"rendering {len(sel)} points (pre-vetted {len(pool)} of {len(cand)}; "
              f"by admitted-view count)")
    else:
        sel = [(pid, T, vis, extra, None) for pid, T, vis, extra in pool[: args.num]]
        print(f"rendering {len(sel)} points (of {len(cand)} with extra coverage)")

    sel_ids = [r[0] for r in sel]
    vis_of = {r[0]: r[2] for r in sel}
    track_of = {r[0]: r[1] for r in sel}

    # Stages 2-4: grow the view set by agreement, re-refining each round.
    n_cur = dict(n_track)
    good = {pid: sorted(track_of[pid]) for pid in sel_ids}
    corr0: dict[int, dict[int, float]] = {}  # round-0 corr to the track template
    thr0: dict[int, float] = {}
    for r in range(args.iters):
        view_idx = [[] for _ in s.cloud.point_ids]
        for pid in sel_ids:
            A = vis_of[pid]
            _, corr, thr, admitted, _ = vet(
                s, pid, n_cur[pid], A, track_of[pid], good[pid], w,
                render_res, ext_of[pid], args.accept_frac, args.floor,
            )
            if r == 0:
                corr0[pid] = {A[k]: float(corr[k]) for k in range(len(A))}
                thr0[pid] = thr
            good[pid] = sorted(track_of[pid] | set(admitted))
            view_idx[cidx[pid]] = [int(i) for i in good[pid]]
        cloud_g = PatchCloud.from_reconstruction(recon, normal="stored", extent_value=EXTENT_FACTOR)
        res_g = cloud_g.refine_normals(
            recon, imgs, point_ids=sel_ids, view_indices=view_idx, **common
        )
        n_cur = {pid: np.asarray(res_g["normal"])[cidx[pid]] for pid in sel_ids}
    n_good = n_cur

    blocks: list[np.ndarray] = []
    for pid in sel_ids:
        A = vis_of[pid]
        T, G = track_of[pid], set(good[pid])
        ext = ext_of[pid]
        accepted = G - T
        rejected = set(A) - G

        # Render and score over the *full* visible set, then optionally subsample
        # which tiles are drawn (decisions/Φ stay over the full set).
        tiles_t, cores_t = render_cores(s, pid, n_track[pid], A, render_res, ext)
        rows_t = znorm_stack(cores_t, w)
        tiles_g, cores_g = render_cores(s, pid, n_good[pid], A, render_res, ext)
        rows_g = znorm_stack(cores_g, w)
        good_pos = [A.index(i) for i in good[pid]]
        wt_g = irls_weights(rows_g[good_pos]) if len(good_pos) >= 2 else np.ones(len(good_pos))
        wmap = {good[pid][k]: float(wt_g[k]) for k in range(len(good_pos))}
        v = len(good_pos)

        if args.max_views and len(A) > args.max_views:
            d = np.linspace(0, len(A) - 1, args.max_views).round().astype(int)
            di = sorted(dict.fromkeys(d.tolist()))
        else:
            di = list(range(len(A)))
        D = [A[k] for k in di]
        cols = [TRACK_COL if i in T else (GOOD_COL if i in G else DROP_COL) for i in D]
        base = float(np.mean([corr0[pid][i] for i in sorted(T)]))

        # Top strip: track normal, each shown view annotated with its correlation
        # to the track template; colour encodes the admit/reject decision.
        c0 = [corr0[pid][i] for i in D]
        bar0 = [float(np.clip(c / base, 0, 1)) if base > 1e-6 else 0.0 for c in c0]
        row_t = strip_image([tiles_t[k] for k in di], D, cols,
                            [f"{c:+.2f}" for c in c0], bar0, render_res, label="track normal")

        # Bottom strip: good-view normal; admitted views carry their inclusion
        # weight, rejected views are greyed out (not in the consensus).
        notes_g = [f"{wmap[i]:.2f}" if i in G else "x" for i in D]
        bars_g = [float(np.clip(wmap[i] * v, 0, 2) / 2) if i in G else None for i in D]
        row_g = strip_image([tiles_g[k] for k in di], D, cols, notes_g, bars_g,
                            render_res, label="good-view normal")

        phi_track = _phi(rows_t[good_pos], irls_weights(rows_t[good_pos]))
        phi_good = _phi(rows_g[good_pos], wt_g)
        dang = _angle(n_track[pid], n_good[pid])

        header = np.full((30, row_t.shape[1], 3), 18, np.uint8)
        shown = f"; showing {len(D)}" if len(D) < len(A) else ""
        msg = (f"pt {pid}  {len(T)} track  +{len(accepted)} admitted  "
               f"-{len(rejected)} rejected  (of {len(A)} visible{shown})   "
               f"thr={thr0[pid]:.2f}   normal moved {dang:.1f} deg   "
               f"consensus Phi over good set: track-n={phi_track:+.3f} -> good-n={phi_good:+.3f}")
        cv2.putText(header, msg, (6, 20), cv2.FONT_HERSHEY_SIMPLEX, 0.44, (180, 220, 180), 1, cv2.LINE_AA)
        blocks.extend([header, row_t, row_g, np.full((6, header.shape[1], 3), 18, np.uint8)])

        print(f"pt {pid:4d}: {len(T):2d} track +{len(accepted):2d} admitted "
              f"-{len(rejected):2d} rejected | thr={thr0[pid]:.2f} | "
              f"normal {dang:5.1f}deg | Phi(good) {phi_track:+.3f}->{phi_good:+.3f}")

    if args.sift_overlap:
        _report_sift_overlap(s, recon, sel_ids, track_of, good)

    width = max(b.shape[1] for b in blocks)
    blocks = [np.pad(b, ((0, 0), (0, width - b.shape[1]), (0, 0)), constant_values=18) for b in blocks]
    legend = legend_bar(width, [
        (TRACK_COL, "track view (seed)"),
        (GOOD_COL, "admitted"),
        (DROP_COL, "rejected (dimmed)"),
        ((255, 255, 255), "validated extent"),
    ])
    montage = np.vstack([legend, *blocks])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(args.output), montage)
    print(f"wrote {args.output} ({montage.shape[1]}x{montage.shape[0]})")


if __name__ == "__main__":
    main()
