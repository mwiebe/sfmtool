# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 10 — what separates real-track clusters from spurious ones?

Run the density-seeded clustering, then split the produced clusters into:
  * track clusters    — contain >=2 members of a single solve point (a real
                        co-observation was captured), and
  * non-track clusters — no two members share a solve point (spurious grouping
                        of background / unrelated features),
and plot a battery of per-cluster descriptor-space statistics for the two
populations, looking for a label-free discriminator.

Statistics per cluster: size, mean/max distance to centroid, spread, mean
pairwise distance, an isolation gap (centroid distance to the first non-member
minus the last member), and the churn after one mean-shift step.

Usage:
    pixi run -e experiments python experiments/exp10_cluster_stats.py seoul_bull_ws
"""

from __future__ import annotations

import argparse
import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from exp05_cluster_match import K, derive_threshold
from seed_cluster import seed_claim_clusters
from sfm_descriptors import load_descriptor_bank

QK = 64


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("workspace")
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    from sfmtool import KdForest

    path = sorted(glob.glob(f"../{args.workspace}/sfmr/*solve*.sfmr"))[0]
    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    img = bank.image_label
    point = bank.point_label
    forest = KdForest(desc, preset="accurate")
    idx, dst = forest.query(desc, k=K)
    idx = idx.astype(np.int64)
    T = derive_threshold(dst[:, 1], "otsu")

    labels = seed_claim_clusters(idx, dst, img, T)
    # group members by cluster id
    valid = np.flatnonzero(labels >= 0)
    order = valid[np.argsort(labels[valid], kind="stable")]
    clab = labels[order]
    bounds = np.r_[0, np.flatnonzero(np.diff(clab)) + 1, len(order)]
    clusters = [order[s:e] for s, e in zip(bounds[:-1], bounds[1:])]

    # centroids + their neighbours (for isolation gap and the mean-shift churn)
    cent = np.zeros((len(clusters), desc.shape[1]), dtype=np.uint8)
    for i, m in enumerate(clusters):
        cent[i] = np.clip(np.round(desc[m].astype(np.float32).mean(0)), 0, 255)
    cidx, cdst = forest.query(cent, k=QK)

    rows = []
    for i, m in enumerate(clusters):
        xc = desc[m].astype(np.float32)
        c = xc.mean(0)
        dctr = np.sqrt(((xc - c) ** 2).sum(1))
        # pairwise mean
        if len(m) <= 64:
            pw = np.sqrt(((xc[:, None] - xc[None]) ** 2).sum(-1))
            mean_pw = pw[np.triu_indices(len(m), 1)].mean() if len(m) > 1 else 0.0
        else:
            mean_pw = 2.0 * dctr.mean()
        # isolation: among the centroid's returned neighbours, last member vs
        # first non-member distance
        mset = set(int(x) for x in m)
        nbr, nbd = cidx[i], cdst[i]
        is_mem = np.array([int(x) in mset for x in nbr])
        max_mem_d = nbd[is_mem].max() if is_mem.any() else 0.0
        nonmem_d = nbd[~is_mem]
        first_non = nonmem_d.min() if nonmem_d.size else np.nan
        gap = first_non - max_mem_d
        # one mean-shift step churn (no claiming)
        mm = nbd <= T
        ci2, cd2 = nbr[mm], nbd[mm]
        if ci2.size:
            o = np.lexsort((cd2, img[ci2]))
            ii = img[ci2][o]
            keep = np.empty(len(o), bool)
            keep[0] = True
            keep[1:] = ii[1:] != ii[:-1]
            new = set(int(x) for x in ci2[o][keep])
        else:
            new = set()
        churn = len(mset ^ new) / max(len(mset | new), 1)
        # label: track if >=2 members share a solve point
        pm = point[m]
        pmi = pm[pm >= 0]
        dom = np.bincount(pmi).max() if pmi.size else 0
        rows.append(
            (dom >= 2, len(m), dctr.mean(), dctr.max(), dctr.std(), mean_pw, gap, churn)
        )

    rows = np.array(rows, dtype=float)
    is_track = rows[:, 0] > 0
    n_tr, n_non = int(is_track.sum()), int((~is_track).sum())
    print(f"{args.workspace}: clusters={len(clusters)} track={n_tr} non-track={n_non}")

    names = [
        "size",
        "mean dist→centroid",
        "max dist→centroid",
        "spread (std)",
        "mean pairwise",
        "isolation gap",
        "mean-shift churn",
    ]
    cols = [1, 2, 3, 4, 5, 6, 7]
    fig, axes = plt.subplots(2, 4, figsize=(18, 8))
    for ax, name, c in zip(axes.ravel(), names, cols):
        tv = rows[is_track, c]
        nv = rows[~is_track, c]
        tv = tv[np.isfinite(tv)]
        nv = nv[np.isfinite(nv)]
        lo = float(min(tv.min(), nv.min()))
        hi = float(np.percentile(np.concatenate([tv, nv]), 99))
        bins = np.linspace(lo, hi, 60)
        ax.hist(tv, bins=bins, density=True, alpha=0.6, label=f"track ({n_tr})")
        ax.hist(nv, bins=bins, density=True, alpha=0.6, label=f"non-track ({n_non})")
        ax.set_title(name)
        ax.legend(fontsize=8)
    axes.ravel()[-1].axis("off")
    fig.suptitle(
        f"{args.workspace}: per-cluster stats, track vs non-track (Otsu T={T:.0f})"
    )
    fig.tight_layout()
    import os

    os.makedirs(args.outdir, exist_ok=True)
    f = f"{args.outdir}/exp10_{args.workspace}_clusterstats.png"
    fig.savefig(f, dpi=110)
    print(f"wrote {f}")

    # quick numeric separation summary (medians)
    print(f"  {'stat':<20} {'track med':>10} {'non med':>10}")
    for name, c in zip(names, cols):
        tv = rows[is_track, c]
        nv = rows[~is_track, c]
        tv = tv[np.isfinite(tv)]
        nv = nv[np.isfinite(nv)]
        print(f"  {name:<20} {np.median(tv):>10.2f} {np.median(nv):>10.2f}")

    # size-stratified: is anything separable at FIXED size (esp. the ambiguous
    # size-2 pairs, where size itself can't help)?
    sz = rows[:, 1]
    for lo, hi, lbl in [(2, 2, "size==2"), (3, 999, "size>=3")]:
        sel = (sz >= lo) & (sz <= hi)
        tr = sel & is_track
        no = sel & ~is_track
        print(f"\n  [{lbl}] track={int(tr.sum())} non-track={int(no.sum())}")
        print(f"    {'stat':<20} {'track med':>10} {'non med':>10}")
        for name, c in zip(names[1:], cols[1:]):
            tv = rows[tr, c]
            nv = rows[no, c]
            tv = tv[np.isfinite(tv)]
            nv = nv[np.isfinite(nv)]
            if tv.size and nv.size:
                print(f"    {name:<20} {np.median(tv):>10.2f} {np.median(nv):>10.2f}")

    # size==2 focused figure (pair distance, isolation gap, churn)
    sel2 = sz == 2
    fig2, ax2 = plt.subplots(1, 3, figsize=(14, 4))
    for ax, name, c in zip(
        ax2, ["mean pairwise", "isolation gap", "mean-shift churn"], [5, 6, 7]
    ):
        tv = rows[sel2 & is_track, c]
        nv = rows[sel2 & ~is_track, c]
        tv = tv[np.isfinite(tv)]
        nv = nv[np.isfinite(nv)]
        hi = float(np.percentile(np.concatenate([tv, nv]), 99))
        lo = float(min(tv.min(), nv.min()))
        bins = np.linspace(lo, hi, 50)
        ax.hist(tv, bins=bins, density=True, alpha=0.6, label=f"track ({tv.size})")
        ax.hist(nv, bins=bins, density=True, alpha=0.6, label=f"non-track ({nv.size})")
        ax.set_title(f"size==2: {name}")
        ax.legend(fontsize=8)
    fig2.suptitle(f"{args.workspace}: size-2 clusters only (the ambiguous pairs)")
    fig2.tight_layout()
    f2 = f"{args.outdir}/exp10_{args.workspace}_size2.png"
    fig2.savefig(f2, dpi=110)
    print(f"wrote {f2}")


if __name__ == "__main__":
    main()
