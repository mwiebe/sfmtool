# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Several candidate visualizations of how the per-point floor relates to tracks.

For each plot idea, emits one figure with a panel per reconstruction (``--all``),
or a single dataset's three figures (``--solve``):

  idea1  recall vs background-admitted as the floor scale (alpha) is swept
  idea2  co-observation vs background distance distributions, with the floor
  idea3  per-descriptor neighbour-distance profiles (co-obs vs background), with
         B and alpha*B marked, across a range of track sizes

Usage:
    pixi run -e experiments python experiments/make_radius_ideas.py --all --outdir /tmp
    pixi run -e experiments python experiments/make_radius_ideas.py \
        --solve ../seattle_backyard_ws/sfmr/*solve*.sfmr --outdir /tmp
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from sfm_descriptors import SOLVES, load_descriptor_bank, resolve_solve
from sfmtool import KdForest

D_RANK = 28
BG_ALPHA = 0.8
KQ = 64

GREEN, RED, DARK, GREY = "#2a9d8f", "#e76f51", "#264653", "#c0c4cc"
SIZES = [2, 4, 8, 16]  # track sizes for the idea3 profile columns


def gather(bank):
    desc = np.ascontiguousarray(bank.descriptors)
    idx, dst = KdForest(desc, preset="accurate").query(desc, k=KQ)
    idx = idx.astype(np.int64)
    B = dst[:, D_RANK].astype(np.float32)
    pt, img = bank.point_label, bank.image_label
    Xf = desc.astype(np.float32)

    # exact co-observation distances per in-track descriptor (directed i->member)
    co_src, co_dist = [], []
    order = np.argsort(pt, kind="stable")
    pts = pt[order]
    bnd = np.r_[0, np.flatnonzero(np.diff(pts)) + 1, len(pts)]
    for s, e in zip(bnd[:-1], bnd[1:]):
        if pts[s] < 0 or e - s < 2:
            continue
        mem = order[s:e]
        D = np.linalg.norm(Xf[mem][:, None] - Xf[mem][None], axis=2)
        iu, ju = np.where(~np.eye(len(mem), dtype=bool))
        co_src.extend(mem[iu].tolist())
        co_dist.extend(D[iu, ju].tolist())

    n = bank.n
    i_rep = np.repeat(np.arange(n), KQ)
    j = idx.ravel()
    dd = dst.ravel().astype(np.float32)
    cross = img[i_rep] != img[j]
    same_track = (pt[i_rep] >= 0) & (pt[i_rep] == pt[j])
    bg = cross & ~same_track & (i_rep != j)
    return dict(
        n=n,
        idx=idx,
        dst=dst,
        B=B,
        co_src=np.array(co_src),
        co_dist=np.array(co_dist, np.float32),
        bg_src=i_rep[bg],
        bg_dist=dd[bg],
        pt=pt,
    )


def draw_idea1(ax, g, title):
    alphas = np.linspace(0.2, 1.4, 49)
    n = g["n"]
    in_track = np.unique(g["co_src"])
    ncov = np.bincount(g["co_src"], minlength=n).astype(float)
    rec, bgc = [], []
    for a in alphas:
        r = a * g["B"]
        per = np.bincount(g["co_src"], g["co_dist"] <= r[g["co_src"]], minlength=n)
        rec.append(np.mean((per / np.maximum(ncov, 1))[in_track]))
        bgc.append((g["bg_dist"] <= r[g["bg_src"]]).sum() / len(in_track))
    ax.plot(alphas, rec, color=GREEN, lw=2)
    ax.set_ylim(0, 1)
    ax.tick_params(axis="y", labelcolor=GREEN)
    ax2 = ax.twinx()
    ax2.plot(alphas, bgc, color=RED, lw=2)
    ax2.tick_params(axis="y", labelcolor=RED)
    ax.axvline(BG_ALPHA, color=DARK, ls="--", lw=1)
    ax.set_title(title, fontsize=9)
    ax.spines[["top"]].set_visible(False)
    ax2.spines[["top"]].set_visible(False)


def draw_idea2(ax, g, title):
    floor = BG_ALPHA * g["B"]
    hi = float(np.percentile(np.r_[g["co_dist"], floor], 99))
    bins = np.linspace(0, hi, 70)
    ax.hist(
        np.clip(g["co_dist"], 0, hi), bins=bins, density=True, color=GREEN, alpha=0.55
    )
    ax.hist(
        np.clip(g["bg_dist"], 0, hi), bins=bins, density=True, color=RED, alpha=0.45
    )
    fm = float(np.median(floor))
    ax.axvline(fm, color=DARK, ls="--", lw=1.5)
    ax.set_xlim(0, hi)
    ax.set_title(f"{title}  (floor ≈ {fm:.0f})", fontsize=9)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)


def draw_idea3_cell(ax, g, ts):
    pt = g["pt"]
    tsize = np.bincount(pt[pt >= 0])
    avail = np.unique(tsize[tsize >= 2])
    ts = int(avail[np.argmin(np.abs(avail - ts))]) if len(avail) else ts
    cand = np.flatnonzero((pt >= 0) & (tsize[np.clip(pt, 0, None)] == ts))
    if len(cand) == 0:
        ax.set_axis_off()
        return None
    i = int(cand[len(cand) // 2])
    ranks = np.arange(1, KQ)
    dd, jj = g["dst"][i, 1:], g["idx"][i, 1:]
    is_co = (g["pt"][jj] >= 0) & (g["pt"][jj] == pt[i])
    ax.scatter(ranks[~is_co], dd[~is_co], s=7, color=GREY)
    ax.scatter(ranks[is_co], dd[is_co], s=11, color=GREEN)
    ax.axhline(g["B"][i], color=DARK, ls=":", lw=1)
    ax.axhline(BG_ALPHA * g["B"][i], color=RED, ls="--", lw=1)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    return ts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--solve", help="single dataset (path/glob)")
    ap.add_argument("--all", action="store_true", help="all four reconstructions")
    ap.add_argument("--outdir", default="/tmp")
    args = ap.parse_args()
    od = Path(args.outdir)
    od.mkdir(parents=True, exist_ok=True)

    datasets = (
        [(lbl, glob) for lbl, glob in SOLVES if lbl != "dino_large"]
        if args.all
        else [("dataset", args.solve)]
    )
    nd = len(datasets)
    rows = (nd + 1) // 2
    f1, ax1 = plt.subplots(rows, 2, figsize=(8.2, 3.0 * rows), squeeze=False)
    f2, ax2 = plt.subplots(rows, 2, figsize=(8.2, 3.0 * rows), squeeze=False)
    f3, ax3 = plt.subplots(
        nd, len(SIZES), figsize=(2.0 * len(SIZES), 2.0 * nd), squeeze=False
    )

    for d, (lbl, glob) in enumerate(datasets):
        print(f"gather {lbl} ...", flush=True)
        g = gather(load_descriptor_bank(resolve_solve(glob)))
        draw_idea1(ax1.ravel()[d], g, lbl)
        draw_idea2(ax2.ravel()[d], g, lbl)
        for c, ts in enumerate(SIZES):
            got = draw_idea3_cell(ax3[d][c], g, ts)
            if d == 0 and got:
                ax3[d][c].set_title(f"track size {got}", fontsize=9)
        ax3[d][0].set_ylabel(lbl, fontsize=8)
        del g

    for d in range(nd, rows * 2):  # hide unused panels
        ax1.ravel()[d].set_axis_off()
        ax2.ravel()[d].set_axis_off()

    f1.suptitle(
        "Recall (green, left) vs background admitted (red, right) vs α; dashed = α0.8",
        fontsize=10,
    )
    f2.suptitle(
        "Co-observation (green) vs background (red) distances; "
        "dashed = median floor α·B",
        fontsize=10,
    )
    f3.suptitle(
        "Neighbour distance vs rank: co-obs (green), background (grey); "
        "α·B dashed, B dotted",
        fontsize=10,
    )
    for f, name in [
        (f1, "all-idea1-alpha-sweep"),
        (f2, "all-idea2-distance-dists"),
        (f3, "all-idea3-profiles"),
    ]:
        f.tight_layout()
        f.savefig(od / f"{name}.png", dpi=120)
        plt.close(f)
    print(f"wrote three figures to {od}")


if __name__ == "__main__":
    main()
