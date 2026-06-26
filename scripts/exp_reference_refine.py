#!/usr/bin/env python3
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment: sharpen the reference patch by congealing the view stack.

`exp_view_localization.py` showed each view's projection carries an incoherent
sub-pixel residual (~0.4 px on clean scenes, several px on the non-convex dino).
This is the follow-up: if we feed each view's **sub-pixel-shifted** render into
the consensus instead of the raw projection, does the reference patch sharpen?

It runs a **congealing** loop (group-wise translation registration):

1. Render every view from the source at its accumulated in-plane offset -- a
   *single* resample (translating the patch center, never re-sampling an
   already-warped tile), so applying a shift never compounds blur.
2. Build the robust consensus.
3. For each view, search the residual shift to the **leave-one-out** consensus
   (a view never aligns to a template its own pixels polluted -- the guard against
   a mean that just fits itself), and add it to that view's offset.
4. Repeat to convergence (mean residual shift -> 0).

Whether the sharpening is *real* is judged by **leave-one-out ZNCC** -- each
view's agreement with the consensus of the *others*. That can only rise if the
views genuinely register better; a mean fitting its own noise would not move it.
Reported alongside the consensus Phi, the consensus-image gradient energy
(sharpness), and the total shift applied. Outputs a before/after consensus-patch
montage.

Usage::

    pixi run python scripts/exp_reference_refine.py RECON.sfmr -o out.png
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
    _template,
    gauss_window,
    irls_weights,
    legend_bar,
    znorm_stack,
)
from exp_view_localization import (  # noqa: E402
    _zncc,
    loo_reference,
    patch_frame,
    render_tile_at_offset,
    search_shift,
    select_good_view_points,
    wpp_to_src,
)

from sfmtool._solve_strips import _SolveStrips  # noqa: E402
from sfmtool._sfmtool import SfmrReconstruction  # noqa: E402


def consensus_image(cores, w):
    """Robust (IRLS) weighted-mean consensus *image* (P, P, C) over the stack."""
    rows = znorm_stack(cores, w)
    wt = irls_weights(rows) if len(cores) >= 2 else np.ones(len(cores))
    raw = np.stack([np.asarray(c, np.float64) for c in cores])
    return (wt[:, None, None, None] * raw).sum(0)


def loo_zncc(cores, w) -> float:
    """Mean over views of ZNCC(view, robust consensus of the *other* views) --
    the honest agreement metric: it rises only if the views truly co-register."""
    rows = znorm_stack(cores, w)
    v = len(cores)
    if v < 2:
        return float("nan")
    vals = []
    for k in range(v):
        idx = [j for j in range(v) if j != k]
        sub = rows[idx]
        wt = (
            irls_weights(sub) if len(idx) >= 2 else np.ones(len(idx)) / max(len(idx), 1)
        )
        xbar = _template(sub, wt)
        nrm = np.sqrt((xbar**2).sum(1, keepdims=True))
        m = xbar / np.where(nrm > 1e-9, nrm, 1.0)
        vals.append(_zncc(rows[k], m))
    return float(np.mean(vals))


def sharpness(img) -> float:
    """Gradient energy of the consensus image -- higher = sharper (less
    registration blur). Same surfel content before/after, so it is comparable."""
    g = img.astype(np.float32)
    g = g.mean(2) if g.ndim == 3 else g
    gx = cv2.Sobel(g, cv2.CV_32F, 1, 0)
    gy = cv2.Sobel(g, cv2.CV_32F, 0, 1)
    return float((gx * gx + gy * gy).mean())


def congeal(s, pid, normal, G, ext, render_res, off, w, *, iters, search, coarse):
    """Congeal the view stack for one point. Returns ``(stats, img0, imgN)``."""
    center, u_ax, v_ax, wpp = patch_frame(s, pid, normal, G, ext, render_res)
    nv = len(G)

    def render(acc):
        tiles = [
            render_tile_at_offset(
                s,
                G[v],
                center,
                normal,
                u_ax,
                v_ax,
                acc[v, 0],
                acc[v, 1],
                wpp,
                ext,
                render_res,
            )
            for v in range(nv)
        ]
        cores = [t[off : off + PATCH, off : off + PATCH] for t in tiles]
        return tiles, cores

    acc = np.zeros((nv, 2))
    tiles, cores = render(acc)
    img0 = consensus_image(cores, w)
    rows0 = znorm_stack(cores, w)
    stat0 = (_phi(rows0, irls_weights(rows0)), loo_zncc(cores, w), sharpness(img0))

    used = 0
    for used in range(1, iters + 1):
        rows = znorm_stack(cores, w)
        deltas = np.zeros((nv, 2))
        for v in range(nv):
            ref = loo_reference(cores, rows, v)
            dx, dy, _, _ = search_shift(tiles[v], ref, off, w, search, coarse)
            deltas[v] = (dx, dy)
        acc = np.clip(acc + deltas, -search, search)  # bound total drift
        tiles, cores = render(acc)
        if float(np.hypot(deltas[:, 0], deltas[:, 1]).mean()) < 0.05:
            break

    imgN = consensus_image(cores, w)
    rowsN = znorm_stack(cores, w)
    statN = (_phi(rowsN, irls_weights(rowsN)), loo_zncc(cores, w), sharpness(imgN))

    mag = np.hypot(acc[:, 0], acc[:, 1])
    src_per_px = wpp_to_src(s, G, center, u_ax, v_ax, wpp)
    stats = {
        "views": nv,
        "iters": used,
        "shift_src": float(np.median(mag)) * src_per_px,
        "phi": (stat0[0], statN[0]),
        "loo": (stat0[1], statN[1]),
        "sharp": (stat0[2], statN[2]),
    }
    return stats, img0, imgN


def _panel(img, label, scale=8):
    """Upscale a consensus patch image to a labelled display tile."""
    p8 = np.clip(img, 0, 255).astype(np.uint8)
    big = cv2.resize(
        p8, (PATCH * scale, PATCH * scale), interpolation=cv2.INTER_NEAREST
    )
    if big.ndim == 2:
        big = cv2.cvtColor(big, cv2.COLOR_GRAY2BGR)
    cv2.rectangle(big, (0, 0), (big.shape[1] - 1, 16), (15, 15, 15), -1)
    cv2.putText(
        big,
        label,
        (4, 12),
        cv2.FONT_HERSHEY_SIMPLEX,
        0.4,
        (230, 230, 230),
        1,
        cv2.LINE_AA,
    )
    return big


def _montage_row(pid, stats, img0, imgN, gutter=170):
    disp = PATCH * 8
    before = _panel(img0, "before (raw projection)")
    after = _panel(imgN, "after (congealed)")
    sep = np.full((disp, 3, 3), 40, np.uint8)
    g = np.full((disp, gutter, 3), 25, np.uint8)
    lines = [
        f"pt {pid}  ({stats['views']} views)",
        f"iters {stats['iters']}",
        f"shift {stats['shift_src']:.2f}px",
        f"LOO {stats['loo'][0]:.3f}->{stats['loo'][1]:.3f}",
        f"Phi {stats['phi'][0]:.3f}->{stats['phi'][1]:.3f}",
        f"sharp x{stats['sharp'][1] / max(stats['sharp'][0], 1e-9):.2f}",
    ]
    for k, t in enumerate(lines):
        cv2.putText(
            g,
            t,
            (6, 26 + 26 * k),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.42,
            (225, 225, 225),
            1,
            cv2.LINE_AA,
        )
    return np.hstack([g, before, sep, after])


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("recon", type=Path, help="solved .sfmr")
    ap.add_argument("-o", "--out", type=Path, default=Path("reference_refine.png"))
    ap.add_argument("--num", type=int, default=8, help="points to analyze")
    ap.add_argument(
        "--montage", type=int, default=4, help="points drawn in the montage"
    )
    ap.add_argument(
        "--context", type=int, default=64, help="context render size (>= PATCH)"
    )
    ap.add_argument("--search", type=int, default=6, help="max total shift (patch px)")
    ap.add_argument("--iters", type=int, default=5, help="congealing iterations")
    ap.add_argument("--no-coarse", action="store_true", help="full-res search only")
    ap.add_argument("--pool", type=int, default=120, help="candidate points pre-vetted")
    ap.add_argument("--accept-frac", type=float, default=0.7)
    ap.add_argument("--floor", type=float, default=0.1)
    args = ap.parse_args()

    recon = SfmrReconstruction.load(str(args.recon))
    s = _SolveStrips(
        recon, recon.workspace_dir, patch=PATCH, extent_factor=EXTENT_FACTOR
    )
    imgs = [s.image(i) for i in range(len(s.names))]
    render_res = max(args.context, PATCH)
    off = (render_res - PATCH) // 2
    w = gauss_window(PATCH)
    common = dict(resolution=PATCH, angular_range_deg=ANG_RANGE, init_steps=INIT_STEPS)
    coarse = not args.no_coarse

    pts = select_good_view_points(
        recon,
        s,
        imgs,
        render_res,
        w,
        common,
        num=args.num,
        pool=args.pool,
        accept_frac=args.accept_frac,
        floor=args.floor,
    )

    print(
        f"\n{'pt':>6} {'views':>5} {'shift':>7} {'LOO ZNCC 0->N':>20} "
        f"{'Phi 0->N':>18} {'sharpness x':>11} {'it':>3}"
    )
    print("-" * 76)

    rows_m, agg = [], {"loo0": [], "looN": [], "dloo": [], "dphi": [], "sharp": []}
    for p in pts:
        stats, img0, imgN = congeal(
            s,
            p["pid"],
            p["normal"],
            p["good"],
            p["ext"],
            render_res,
            off,
            w,
            iters=args.iters,
            search=args.search,
            coarse=coarse,
        )
        l0, lN = stats["loo"]
        ph0, phN = stats["phi"]
        sr = stats["sharp"][1] / max(stats["sharp"][0], 1e-9)
        agg["loo0"].append(l0)
        agg["looN"].append(lN)
        agg["dloo"].append(lN - l0)
        agg["dphi"].append(phN - ph0)
        agg["sharp"].append(sr)
        print(
            f"{p['pid']:>6} {stats['views']:>5} {stats['shift_src']:>5.2f}px  "
            f"{l0:>+7.3f}->{lN:<+7.3f}  {ph0:>+6.3f}->{phN:<+6.3f}  "
            f"x{sr:>6.2f}  {stats['iters']:>3}"
        )
        if len(rows_m) < args.montage:
            rows_m.append(_montage_row(p["pid"], stats, img0, imgN))

    print("-" * 76)
    print(
        f"median over {len(pts)} pts: LOO ZNCC {np.median(agg['loo0']):+.3f}->"
        f"{np.median(agg['looN']):+.3f} (dLOO {np.median(agg['dloo']):+.3f}), "
        f"dPhi {np.median(agg['dphi']):+.3f}, sharpness x{np.median(agg['sharp']):.2f}"
    )

    if rows_m:
        width = max(r.shape[1] for r in rows_m)
        rows_m = [
            np.hstack([r, np.full((r.shape[0], width - r.shape[1], 3), 25, np.uint8)])
            if r.shape[1] < width
            else r
            for r in rows_m
        ]
        legend = legend_bar(
            width, [((230, 230, 230), "consensus patch: before vs after congealing")]
        )
        out = np.vstack([legend, *rows_m])
        args.out.parent.mkdir(parents=True, exist_ok=True)
        cv2.imwrite(str(args.out), out)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
