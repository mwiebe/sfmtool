#!/usr/bin/env python3
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Assess whether good-view refinement makes the normals more robust/correct.

Two complementary read-outs, both comparing **track-only** vs. **good-view**
refinement on the *same* points:

1. **Quantitative robustness.**
   * *Confidence* — the refiner's Φ-peakedness (Hessian of the photoconsistency
     at the optimum); a sharper peak means a better-determined normal.
   * *Leave-one-out stability* — re-refine the normal dropping each view in turn
     and report the angular spread of the results. A 3-view track has no
     redundancy (drop one ⇒ below `min_views`), so it cannot even be evaluated;
     the good-view set, backed by many vetted views, should be stable to a degree
     or two. This is the crux of "robustness".

2. **Visual check.** Writes two reconstructions —
   ``<stem>.track.sfmr`` and ``<stem>.good.sfmr`` — each carrying the refined
   patch cloud and RGBA patch bitmaps, ready for::

       sfm render-patches <stem>.good.sfmr -o renders/ --mode texture --opaque

   A correct surfel's texture *continues* the underlying image; mis-oriented
   normals smear it. Compare the two overlays (and ``--mode normal`` for the
   normal field) to see the improvement directly on the images.

Usage::

    pixi run python scripts/exp_goodview_render_patches.py RECON.sfmr -o out/stem
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
from sfmtool._sfmtool import PatchCloud, SfmrReconstruction  # noqa: E402


def _angle_deg(a, b) -> float:
    a = a / (np.linalg.norm(a) + 1e-12)
    b = b / (np.linalg.norm(b) + 1e-12)
    return float(np.degrees(np.arccos(np.clip(abs(np.dot(a, b)), -1.0, 1.0))))


def _refine(recon, imgs, point_ids, *, view_indices=None, sampler="bilinear",
            cache="fronto", extent=EXTENT_FACTOR, bitmaps=True):
    cloud = PatchCloud.from_reconstruction(recon, normal="stored", extent_value=extent)
    res = cloud.refine_normals(
        recon, imgs, resolution=PATCH, angular_range_deg=ANG_RANGE, init_steps=INIT_STEPS,
        point_ids=point_ids, view_indices=view_indices, compute_confidence=True,
        render_bitmaps=bitmaps, sampler=sampler, cache=cache,
    )
    return cloud, res


def _loo_spread(recon, imgs, npatch, cidx, pid, views, *, sampler="bilinear",
                cache="fronto", extent=EXTENT_FACTOR) -> float:
    """Angular std (deg) of the normal when each view is dropped in turn."""
    if len(views) < 4:  # need >= min_views(3) after dropping one
        return float("nan")
    normals = []
    for drop in range(len(views)):
        vi = [[] for _ in range(npatch)]
        vi[cidx[pid]] = [int(v) for k, v in enumerate(views) if k != drop]
        cloud = PatchCloud.from_reconstruction(recon, normal="stored", extent_value=extent)
        r = cloud.refine_normals(
            recon, imgs, resolution=PATCH, angular_range_deg=ANG_RANGE,
            init_steps=INIT_STEPS, point_ids=[pid], view_indices=vi,
            sampler=sampler, cache=cache,
        )
        normals.append(np.asarray(r["normal"])[cidx[pid]])
    mean = np.mean(normals, axis=0)
    mean /= np.linalg.norm(mean) + 1e-12
    return float(np.sqrt(np.mean([_angle_deg(n, mean) ** 2 for n in normals])))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("recon", type=Path)
    ap.add_argument("-o", "--out", type=Path, required=True, help="output path stem")
    ap.add_argument("--iters", type=int, default=2)
    ap.add_argument("--accept-frac", type=float, default=0.7)
    ap.add_argument("--floor", type=float, default=0.1)
    ap.add_argument("--loo", type=int, default=12, help="points to leave-one-out (0 = skip)")
    ap.add_argument("--sampler", choices=("bilinear", "anisotropic"), default="bilinear",
                    help="pyramid sampler; anisotropic de-aliases oblique views (slower)")
    ap.add_argument("--cache", choices=("fronto", "off"), default="fronto",
                    help="'off' re-renders every candidate from source (exact, slower)")
    ap.add_argument("--extent", type=float, default=EXTENT_FACTOR,
                    help="patch half-extent as a multiple of feature size (default 5)")
    args = ap.parse_args()
    fid = dict(sampler=args.sampler, cache=args.cache, extent=args.extent)
    print(f"fidelity: sampler={args.sampler}  cache={args.cache}  extent={args.extent}")

    recon = SfmrReconstruction.load(str(args.recon))
    s = _SolveStrips(recon, recon.workspace_dir, patch=PATCH, extent_factor=args.extent)
    imgs = [s.image(i) for i in range(len(s.names))]
    w = gauss_window(PATCH)
    npatch = len(s.cloud)
    cidx = {int(p): k for k, p in enumerate(s.cloud.point_ids)}
    finite = np.abs(np.asarray(recon.positions_xyzw)[:, 3]) > 1e-9

    # Eligible points: finite, track >= 3, and with extra geometrically-visible
    # views to vet. Refine exactly these both ways so the comparison is the normal.
    elig, vis_of, track_of = [], {}, {}
    for pid in s.cloud.point_ids:
        pid = int(pid)
        if not finite[pid]:
            continue
        tracked = sorted(set(s.obs.get(pid, [])))
        if len(tracked) < 3:
            continue
        vis = geometric_views(s, pid)
        if set(vis) - set(tracked):
            elig.append(pid)
            vis_of[pid], track_of[pid] = vis, set(tracked)
    if not elig:
        raise SystemExit("no eligible points")
    print(f"{len(elig)} eligible points (finite, track>=3, with extra visible views)")

    # Track-only refinement of the whole eligible set.
    cloud_t, res_t = _refine(recon, imgs, elig, **fid)
    n_track = {pid: np.asarray(res_t["normal"])[cidx[pid]] for pid in elig}

    # Grow each point's view set by photometric agreement (vet at patch res), then
    # refine the whole eligible set over the good views.
    n_cur = dict(n_track)
    good = {pid: sorted(track_of[pid]) for pid in elig}
    for _ in range(args.iters):
        for pid in elig:
            A = vis_of[pid]
            ext = s._half(pid)
            _, _, _, admitted, _ = vet(
                s, pid, n_cur[pid], A, track_of[pid], good[pid], w,
                PATCH, ext, args.accept_frac, args.floor,
            )
            good[pid] = sorted(track_of[pid] | set(admitted))
        view_idx = [[] for _ in range(npatch)]
        for pid in elig:
            view_idx[cidx[pid]] = [int(i) for i in good[pid]]
        cloud_g, res_g = _refine(recon, imgs, elig, view_indices=view_idx, **fid)
        n_cur = {pid: np.asarray(res_g["normal"])[cidx[pid]] for pid in elig}

    # --- Quantitative robustness ---------------------------------------------
    ct = np.asarray(res_t["confidence"])[[cidx[p] for p in elig]]
    cg = np.asarray(res_g["confidence"])[[cidx[p] for p in elig]]
    ct, cg = ct[np.isfinite(ct)], cg[np.isfinite(cg)]
    dmove = np.array([_angle_deg(n_track[p], n_cur[p]) for p in elig])
    nv = np.array([len(good[p]) for p in elig])
    print("\n=== robustness ===")
    print(f"  views per point:   track 3  ->  good-view {nv.mean():.1f} (median {int(np.median(nv))})")
    print(f"  Phi-peakedness confidence: track {np.nanmean(ct):.3f}  ->  good {np.nanmean(cg):.3f} "
          f"(+{100 * (np.nanmean(cg) / max(np.nanmean(ct), 1e-6) - 1):.0f}%)")
    print(f"  normal moved track->good: median {np.median(dmove):.1f} deg, "
          f"90th pct {np.percentile(dmove, 90):.1f} deg")

    if args.loo:
        rng = np.random.default_rng(0)
        sample = rng.choice(elig, size=min(args.loo, len(elig)), replace=False)
        loo = [_loo_spread(recon, imgs, npatch, cidx, int(p), good[int(p)], **fid)
               for p in sample]
        loo = np.array([x for x in loo if np.isfinite(x)])
        if len(loo):
            print(f"  leave-one-out normal spread (good-view, n={len(loo)}): "
                  f"median {np.median(loo):.2f} deg, max {loo.max():.2f} deg")
        print("  (track has only 3 views: dropping one falls below min_views, so its "
              "normal has no redundancy to test — every view is load-bearing.)")

    # --- Visual: write both reconstructions for render-patches ----------------
    out_t = args.out.with_suffix(".track.sfmr")
    out_g = args.out.with_suffix(".good.sfmr")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    recon.clone_with_changes(patches=cloud_t, patch_bitmaps=res_t["bitmaps"]).save(str(out_t))
    recon.clone_with_changes(patches=cloud_g, patch_bitmaps=res_g["bitmaps"]).save(str(out_g))
    print(f"\nwrote {out_t}\n      {out_g}")
    print("render with, e.g.:")
    print(f"  sfm render-patches {out_g} -o renders_good/ --mode texture --opaque")
    print(f"  sfm render-patches {out_t} -o renders_track/ --mode texture --opaque")


if __name__ == "__main__":
    main()
