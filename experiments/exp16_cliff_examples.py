# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 16 — walk through individual cliff decisions.

Picks a few representative in-track descriptors and plots, for each, its 16
nearest-neighbour distances by rank — coloured by whether the neighbour is a true
co-observation (same 3-D point, different image), cross-image background, or a
same-image feature. Overlaid: `d2` (the reference), and where the 7-nearest cliff
cuts (the largest jump `d[i+1]/d[i]` among the first 7; everything before the jump
is kept). This shows the actual distances and ratios that drive each decision.

Usage:
    pixi run -e experiments python experiments/exp16_cliff_examples.py [dataset_ws]
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np

from exp03_radius_clusters import K
from sfm_descriptors import load_descriptor_bank


def cliff_of(row_dst: np.ndarray):
    """Largest jump among the first 7 neighbours. Returns (cut_rank, outer_rank,
    ratio): cut_rank = last kept neighbour (1-based), outer_rank = first excluded,
    ratio = the jump size d[outer]/d[cut]."""
    d = row_dst[1:8]  # ranks 1..7
    r = d[1:] / np.maximum(d[:-1], 1e-6)
    i = int(np.argmax(r))  # jump between rank i+1 and i+2
    return i + 1, i + 2, float(r[i])


def pick_examples(dst, idx, img, point, pos):
    """Choose 6 illustrative in-track descriptors covering distinct behaviours."""
    rows = np.arange(len(dst))
    nb_idx = idx[:, 1:]
    cross = img[nb_idx] != img[:, None]
    coobs = cross & (point[nb_idx] == point[:, None]) & pos[:, None]
    bg = cross & ~((point[nb_idx] == point[:, None]) & pos[:, None])

    g7 = (dst[:, 2:8] / np.maximum(dst[:, 1:7], 1e-6)).argmax(1)
    outer7 = dst[rows, 2 + g7]
    within = dst[:, 1:] < outer7[:, None]
    co_cnt = coobs.sum(1)
    co_lost = (coobs & ~within).sum(1)
    bg_in = (bg & within).sum(1)

    elig = pos & (co_cnt >= 2)

    def first(mask, key=None, largest=True):
        c = np.flatnonzero(mask & elig)
        if c.size == 0:
            return None
        if key is None:
            return int(c[0])
        v = key[c]
        return int(c[v.argmax() if largest else v.argmin()])

    picks = []
    seen = set()

    def add(s, label):
        if s is not None and s not in seen:
            seen.add(s)
            picks.append((s, label))

    add(first((co_lost == 0) & (bg_in == 0), co_cnt), "clean cut, keeps all co-obs")
    add(first(g7 == 0, co_cnt), "biggest jump is d1→d2 (cuts at rank 1)")
    add(first(co_lost >= 2, co_lost), "intra-track gap (co-obs lost past the cut)")
    add(first(bg_in >= 1, bg_in), "background admitted inside the cut")
    # Fill remaining slots with median-multiplicity in-track points.
    extra = np.flatnonzero(elig)
    extra = extra[np.argsort(co_cnt[extra])][len(extra) // 2 :]
    for s in extra:
        if len(picks) >= 6:
            break
        add(int(s), "typical in-track point")
    return picks


def plot(dst, idx, img, point, pos, picks, name, outpath):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    ranks = np.arange(1, 17)
    for ax, (s, label) in zip(axes.ravel(), picks):
        dk = dst[s, 1:17]
        nb = idx[s, 1:17]
        same_pt = (point[nb] == point[s]) & (point[s] >= 0)
        cross = img[nb] != img[s]
        kind = np.where(same_pt & cross, "co", np.where(cross, "bg", "same"))
        colors = {"co": "C2", "bg": "C1", "same": "0.6"}
        cut_rank, outer_rank, ratio = cliff_of(dst[s])
        d2 = dst[s, 2]
        outer = dst[s, outer_rank]

        ax.plot(ranks, dk, "-", color="0.8", zorder=1)
        for k in ("co", "bg", "same"):
            m = kind == k
            lbl = {"co": "co-observation", "bg": "background", "same": "same image"}[k]
            ax.scatter(ranks[m], dk[m], c=colors[k], s=45, zorder=3, label=lbl)
        ax.axhline(d2, color="C0", ls=":", zorder=2, label="d2")
        ax.axhline(outer, color="k", ls="--", zorder=2, label="cliff cut")
        ax.axhspan(0, outer, color="C2", alpha=0.06, zorder=0)
        # mark the jump that defines the cliff
        ax.annotate(
            f"×{ratio:.2f}",
            xy=(cut_rank + 0.5, (dst[s, cut_rank] + outer) / 2),
            fontsize=9,
            color="k",
            ha="left",
            va="center",
        )
        kept = int((dk < outer).sum())
        ax.set_title(f"{label}\n(cut after rank {cut_rank}; {kept} kept)", fontsize=9)
        ax.set_xlabel("neighbour rank")
        ax.set_ylabel("distance")
        ax.set_xticks(ranks[::2])
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=7, loc="upper left")
    fig.suptitle(f"{name}: how the cliff decides, per descriptor", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, dpi=130)
    print(f"wrote {outpath}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("workspace", nargs="?", default="seoul_bull_ws")
    ap.add_argument("--preset", default="accurate")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    from sfmtool import KdForest

    path = sorted(glob.glob(f"../{args.workspace}/sfmr/*solve*.sfmr"))[0]
    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    idx, dst = KdForest(desc, preset=args.preset).query(desc, k=K)
    idx = idx.astype(np.int64)
    img = bank.image_label.astype(np.int64)
    point = bank.point_label
    pos = point >= 0

    picks = pick_examples(dst, idx, img, point, pos)
    name = args.workspace.replace("_ws", "")
    for s, label in picks:
        cut_rank, outer_rank, ratio = cliff_of(dst[s])
        print(
            f"  #{s}: {label} | d1={dst[s, 1]:.0f} d2={dst[s, 2]:.0f} "
            f"cut@rank{cut_rank} jump×{ratio:.2f} outer={dst[s, outer_rank]:.0f}"
        )
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    plot(
        dst,
        idx,
        img,
        point,
        pos,
        picks,
        name,
        str(out / f"exp16_cliff_examples_{name}.png"),
    )


if __name__ == "__main__":
    main()
