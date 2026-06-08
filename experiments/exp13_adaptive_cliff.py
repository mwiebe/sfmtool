# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 13 — per-point adaptive cliff radius (vs one global T).

exp12 showed the per-point neighbour cliff finds each point's local cluster edge,
but collapsing it to one global radius just re-derives Otsu. Here we keep the
radius *per point*: each seed gathers the neighbours inside its *own* cliff
distance T_i, so dense and sparse regions of descriptor space get tighter/looser
radii automatically. Only points with a real cliff (jump >= 1.3) may seed.

We (1) characterise the raw adaptive-radius data — how T_i is distributed, how it
splits in-track vs background, and how much it drifts between images — and
(2) score adaptive-cliff clusters against the solve's tracks, head-to-head with
the global Otsu(d1) clusterer. A global post-process over T_i is left for later;
this is the raw evaluation.

Usage:
    pixi run -e experiments python experiments/exp13_adaptive_cliff.py
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np

from exp03_radius_clusters import K, score
from exp05_cluster_match import derive_threshold
from exp12_cliff_threshold import cliff
from seed_cluster import seed_claim_clusters
from sfm_descriptors import load_descriptor_bank

DATASETS = [
    "seoul_bull_ws",
    "seattle_backyard_ws",
    "kerry_park_ws",
    "dino_dog_toy_ws",
]
JUMP = 1.3  # cliff-strength gate: keep points whose cliff jump >= this


def one_per_image(members, mdist, img):
    """Resolve members to one feature per image, nearest first."""
    if members.size == 0:
        return members
    o = np.lexsort((mdist, img[members]))
    mi = img[members][o]
    keep = np.empty(len(o), dtype=bool)
    keep[0] = True
    keep[1:] = mi[1:] != mi[:-1]
    return members[o][keep]


def adaptive_clusters(idx, dst, img, ti, can_seed, min_size=2):
    """Density-seeded clustering where seed s gathers within its own radius t_i."""
    n = idx.shape[0]
    rng = np.arange(n)
    valid = (dst <= ti[:, None]) & (idx != rng[:, None]) & (img[idx] != img[:, None])
    order = np.argsort(-valid.sum(1), kind="stable")
    claimed = np.zeros(n, dtype=bool)
    labels = np.full(n, -1, dtype=np.int64)
    cid = 0
    for s in order:
        if claimed[s] or not can_seed[s]:
            continue
        m = valid[s] & ~claimed[idx[s]]
        members = np.concatenate(([s], idx[s][m]))
        mdist = np.concatenate(([0.0], dst[s][m]))
        members = one_per_image(members, mdist, img)
        if members.size >= min_size:
            claimed[members] = True
            labels[members] = cid
            cid += 1
        else:
            claimed[s] = True
    return labels


def score_labels(bank, labels, T):
    """Score a cluster labelling (-1 = unclustered) with exp03's scorer."""
    n = bank.n
    comp = labels.copy()
    unl = np.flatnonzero(labels < 0)
    nxt = int(labels.max()) + 1 if (labels >= 0).any() else 0
    comp[unl] = np.arange(len(unl), dtype=np.int64) + nxt
    sizes = np.bincount(comp, minlength=n)
    return score(bank, comp, sizes, np.ones(n, dtype=bool), T)


def run_one(ws: str, preset: str) -> dict:
    from sfmtool import KdForest

    path = sorted(glob.glob(f"../{ws}/sfmr/*solve*.sfmr"))[0]
    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    idx, dst = KdForest(desc, preset=preset).query(desc, k=K)
    idx = idx.astype(np.int64)
    img = bank.image_label.astype(np.int64)
    pos = bank.point_label >= 0
    d1 = dst[:, 1]

    # Per-point cliff radius from the first 8 neighbours; strength gate.
    cdist, crank, jump = cliff(dst[:, 1:8])
    strong = jump >= JUMP
    T_otsu = derive_threshold(d1, "otsu")

    # Two clusterings: global Otsu(d1) vs per-point adaptive cliff.
    lab_g = seed_claim_clusters(idx, dst, img, T_otsu)
    lab_a = adaptive_clusters(idx, dst, img, cdist, strong)
    sc_g = score_labels(bank, lab_g, T_otsu)
    sc_a = score_labels(bank, lab_a, float(np.median(cdist[strong])))

    name = ws.replace("_ws", "")
    return dict(
        name=name,
        n=bank.n,
        pos=pos,
        strong=strong,
        cdist=cdist,
        jump=jump,
        image=img,
        T_otsu=float(T_otsu),
        sc_g=sc_g,
        sc_a=sc_a,
    )


def print_table(r: dict) -> None:
    p, s = r["pos"], r["strong"]
    print(
        f"\n{r['name']}: n={r['n']}  in-track={p.mean():.1%}  "
        f"strong-cliff={s.mean():.1%}  gate-recall(in-trk strong)={p[s].mean():.1%}"
    )
    print(
        f"  global Otsu T={r['T_otsu']:.0f}  "
        f"adaptive median T_i(strong)={np.median(r['cdist'][s]):.0f}  "
        f"IQR=[{np.percentile(r['cdist'][s], 25):.0f},"
        f"{np.percentile(r['cdist'][s], 75):.0f}]"
    )
    hdr = (
        f"  {'method':<12} {'clusters':>8} {'purity':>7} {'falseM':>7} "
        f"{'recov':>6} {'frags':>6} {'cover':>6}"
    )
    print(hdr)
    for lbl, sc in [("global T", r["sc_g"]), ("adaptive", r["sc_a"])]:
        print(
            f"  {lbl:<12} {sc['components']:>8} {sc['dom_cov']:>7.3f} "
            f"{sc['false_merge_frac']:>7.1%} {sc['recovery']:>6.3f} "
            f"{sc['fragments_per_track']:>6.2f} {sc['in_track_coverage']:>6.1%}"
        )


def plot_radius_hist(results, outpath):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, r in zip(axes.ravel(), results):
        s = r["strong"]
        ci, bg = r["cdist"][s & r["pos"]], r["cdist"][s & ~r["pos"]]
        hi = np.percentile(r["cdist"][s], 99)
        bins = np.linspace(0, hi, 60)
        ax.hist(ci, bins=bins, density=True, alpha=0.6, label="in-track")
        ax.hist(bg, bins=bins, density=True, alpha=0.6, label="background")
        ax.axvline(r["T_otsu"], color="k", ls="--", label=f"global T={r['T_otsu']:.0f}")
        ax.set_title(r["name"])
        ax.set_xlabel("per-point cliff radius T_i (strong cliffs)")
        ax.set_ylabel("density")
        ax.grid(alpha=0.3)
    axes[0, 0].legend()
    fig.suptitle("Per-point adaptive cliff radius vs the single global T", y=0.997)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, dpi=130)
    print(f"\nwrote {outpath}")


def plot_per_image(results, outpath):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, r in zip(axes.ravel(), results):
        m = r["strong"] & r["pos"]
        imgs = r["image"][m]
        ti = r["cdist"][m]
        uimg = np.unique(imgs)
        med = np.array([np.median(ti[imgs == u]) for u in uimg])
        q1 = np.array([np.percentile(ti[imgs == u], 25) for u in uimg])
        q3 = np.array([np.percentile(ti[imgs == u], 75) for u in uimg])
        x = np.arange(len(uimg))
        ax.fill_between(x, q1, q3, alpha=0.25, color="C0", label="IQR")
        ax.plot(x, med, "o-", ms=3, color="C0", label="median T_i")
        ax.axhline(r["T_otsu"], color="k", ls="--", label=f"global T={r['T_otsu']:.0f}")
        ax.set_title(r["name"])
        ax.set_xlabel("image index")
        ax.set_ylabel("in-track cliff radius T_i")
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(
        "Adaptive cliff radius drifts between images (one global T can't)", y=0.997
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, dpi=130)
    print(f"wrote {outpath}")


def plot_metrics(results, outpath):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    keys = [
        ("dom_cov", "purity"),
        ("recovery", "recov"),
        ("in_track_coverage", "cover"),
        ("fragments_per_track", "frags/4"),
    ]
    for ax, r in zip(axes.ravel(), results):
        gv = [r["sc_g"][k] / (4 if k == "fragments_per_track" else 1) for k, _ in keys]
        av = [r["sc_a"][k] / (4 if k == "fragments_per_track" else 1) for k, _ in keys]
        x = np.arange(len(keys))
        ax.bar(x - 0.2, gv, 0.4, label="global T")
        ax.bar(x + 0.2, av, 0.4, label="adaptive")
        ax.set_xticks(x)
        ax.set_xticklabels([lbl for _, lbl in keys])
        ax.set_title(r["name"])
        ax.set_ylim(0, 1.05)
        ax.grid(alpha=0.3, axis="y")
    axes[0, 0].legend()
    fig.suptitle("Cluster quality: global T vs per-point adaptive cliff", y=0.997)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, dpi=130)
    print(f"wrote {outpath}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="accurate")
    ap.add_argument("--datasets", nargs="*", default=DATASETS)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    print(
        "Per-point adaptive cliff radius: raw data + clusters vs global Otsu(d1).\n"
        "purity = dominant-track coverage; falseM = clusters mixing tracks;\n"
        "recov = track completeness; frags = clusters per track; cover = in-track\n"
        "descriptors landed in a multi-cluster."
    )
    results = [run_one(ws, args.preset) for ws in args.datasets]
    for r in results:
        print_table(r)
    if len(results) == len(DATASETS):
        out = Path(args.outdir)
        out.mkdir(parents=True, exist_ok=True)
        plot_radius_hist(results, str(out / "exp13_radius_hist.png"))
        plot_per_image(results, str(out / "exp13_per_image.png"))
        plot_metrics(results, str(out / "exp13_metrics.png"))


if __name__ == "__main__":
    main()
