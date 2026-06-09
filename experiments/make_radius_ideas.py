# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Several candidate visualizations of how the per-point floor relates to tracks.

For one dataset, emits a few different plot ideas (separate PNGs) so we can pick
the clearest:

  idea1  recall vs background-admitted as the floor scale (alpha) is swept
  idea2  co-observation vs background distance distributions, with the floor
  idea3  per-descriptor neighbour-distance profiles (co-obs vs background) for a
         handful of example descriptors, with B and alpha*B marked

Usage:
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

from sfm_descriptors import load_descriptor_bank, resolve_solve
from sfmtool import KdForest

D_RANK = 28
BG_ALPHA = 0.8
KQ = 64

GREEN, RED, DARK = "#2a9d8f", "#e76f51", "#264653"


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
        for a in range(len(mem)):
            for b in range(len(mem)):
                if a != b:
                    co_src.append(int(mem[a]))
                    co_dist.append(float(D[a, b]))
    co_src = np.array(co_src)
    co_dist = np.array(co_dist, np.float32)

    # background neighbours (cross-image, not same track) from the kNN table
    n = bank.n
    i_rep = np.repeat(np.arange(n), KQ)
    j = idx.ravel()
    dd = dst.ravel().astype(np.float32)
    cross = img[i_rep] != img[j]
    same_track = (pt[i_rep] >= 0) & (pt[i_rep] == pt[j])
    bg = cross & ~same_track & (i_rep != j)
    return dict(
        idx=idx,
        dst=dst,
        B=B,
        co_src=co_src,
        co_dist=co_dist,
        bg_src=i_rep[bg],
        bg_dist=dd[bg],
        pt=pt,
        img=img,
    )


def idea1(g, n, out):
    """Recall and background admitted vs the floor scale alpha."""
    alphas = np.linspace(0.2, 1.4, 49)
    in_track = np.unique(g["co_src"])
    ncov = np.bincount(g["co_src"], minlength=n).astype(float)
    nseeds = len(in_track)
    rec, bgc = [], []
    for a in alphas:
        r = a * g["B"]
        cap = g["co_dist"] <= r[g["co_src"]]
        per = np.bincount(g["co_src"], cap, minlength=n)
        rec.append(np.mean((per / np.maximum(ncov, 1))[in_track]))
        capb = g["bg_dist"] <= r[g["bg_src"]]
        bgc.append(capb.sum() / nseeds)
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.plot(alphas, rec, color=GREEN, lw=2, label="co-obs recall")
    ax.set_xlabel("floor scale α  (radius = α·B)")
    ax.set_ylabel("co-observation recall", color=GREEN)
    ax.set_ylim(0, 1)
    ax2 = ax.twinx()
    ax2.plot(alphas, bgc, color=RED, lw=2, label="background admitted")
    ax2.set_ylabel("background neighbours admitted / seed", color=RED)
    ax.axvline(BG_ALPHA, color=DARK, ls="--", lw=1)
    ax.text(BG_ALPHA, 0.02, " α=0.8", color=DARK, fontsize=8)
    for s in ("top",):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def idea2(g, out):
    """Co-observation vs background distance distributions, with the floor."""
    floor = BG_ALPHA * g["B"]
    hi = float(np.percentile(np.r_[g["co_dist"], floor], 99))
    bins = np.linspace(0, hi, 70)
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax.hist(
        np.clip(g["co_dist"], 0, hi),
        bins=bins,
        density=True,
        color=GREEN,
        alpha=0.55,
        label="co-observation distances",
    )
    ax.hist(
        np.clip(g["bg_dist"], 0, hi),
        bins=bins,
        density=True,
        color=RED,
        alpha=0.45,
        label="background distances",
    )
    fm = float(np.median(floor))
    ax.axvline(fm, color=DARK, ls="--", lw=1.5, label=f"median floor α·B ≈ {fm:.0f}")
    ax.set_xlabel("descriptor distance")
    ax.set_ylabel("density")
    ax.set_xlim(0, hi)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def idea3(g, bank, out):
    """Per-descriptor neighbour-distance profiles for example descriptors."""
    pt = g["pt"]
    tsize = np.bincount(pt[pt >= 0])
    targets = [2, 3, 5, 8, 13, 20]
    fig, axes = plt.subplots(2, 3, figsize=(7.6, 4.2), sharex=True)
    for ax, ts in zip(axes.ravel(), targets):
        cand = np.flatnonzero((pt >= 0) & (tsize[np.clip(pt, 0, None)] == ts))
        if len(cand) == 0:
            sizes_av = tsize[tsize >= 2]
            ts = int(sizes_av[np.argmin(np.abs(sizes_av - ts))])
            cand = np.flatnonzero((pt >= 0) & (tsize[np.clip(pt, 0, None)] == ts))
        i = int(cand[len(cand) // 2])
        ranks = np.arange(1, KQ)
        dd = g["dst"][i, 1:]
        jj = g["idx"][i, 1:]
        is_co = (g["pt"][jj] >= 0) & (g["pt"][jj] == pt[i])
        ax.scatter(ranks[~is_co], dd[~is_co], s=10, color="#c0c4cc", label="background")
        ax.scatter(ranks[is_co], dd[is_co], s=14, color=GREEN, label="co-obs")
        B = g["B"][i]
        ax.axhline(B, color=DARK, ls=":", lw=1)
        ax.axhline(BG_ALPHA * B, color=RED, ls="--", lw=1)
        ax.set_title(f"track size {ts}", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(frameon=False, fontsize=7, loc="upper left")
    fig.supxlabel("neighbour rank", fontsize=9)
    fig.supylabel("distance", fontsize=9)
    fig.suptitle("α·B (red dashed) vs B at rank 28 (dotted); co-obs green", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=120)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--solve", required=True)
    ap.add_argument("--outdir", default="/tmp")
    args = ap.parse_args()
    bank = load_descriptor_bank(resolve_solve(args.solve))
    g = gather(bank)
    od = Path(args.outdir)
    idea1(g, bank.n, od / "idea1-alpha-sweep.png")
    idea2(g, od / "idea2-distance-dists.png")
    idea3(g, bank, od / "idea3-profiles.png")
    print(f"wrote idea1/2/3 to {od}")


if __name__ == "__main__":
    main()
