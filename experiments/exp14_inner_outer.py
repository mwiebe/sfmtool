# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 14 — inner radius vs the first non-cluster point.

The adaptive cliff cuts each point's cluster at the inner radius d_r (the last
near member before the biggest relative jump). The *first non-cluster point* is
the very next neighbour, d_{r+1} — there is nothing between them, so d_{r+1} is
exactly what a slightly looser radius would admit next. Comparing the two
histograms answers the question exp13 raised: is the cliff cutting at background
(d_{r+1} is far / wrong-track => good cut) or is it amputating real members
(d_{r+1} is a same-track cross-image co-observation => cut too early)?

For in-track points with a strong cliff we plot the inner radius d_r against the
first-excluded distance d_{r+1}, split by whether that first-excluded neighbour
is actually a same-track co-observation or background, with the global Otsu(d1)
radius marked.

Usage:
    pixi run -e experiments python experiments/exp14_inner_outer.py
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np

from exp03_radius_clusters import K
from exp05_cluster_match import derive_threshold
from exp11_ratio_sweep import youden
from sfm_descriptors import load_descriptor_bank

DATASETS = [
    "seoul_bull_ws",
    "seattle_backyard_ws",
    "kerry_park_ws",
    "dino_dog_toy_ws",
]
JUMP = 1.3
M = 8  # neighbours scanned for the cliff (cols 1..7 are the 7 nearest others)


def run_one(ws: str, preset: str) -> dict:
    from sfmtool import KdForest

    path = sorted(glob.glob(f"../{ws}/sfmr/*solve*.sfmr"))[0]
    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    idx, dst = KdForest(desc, preset=preset).query(desc, k=K)
    idx = idx.astype(np.int64)
    img = bank.image_label.astype(np.int64)
    point = bank.point_label
    pos = point >= 0
    rows = np.arange(bank.n)

    # Cliff over the first M neighbours: biggest relative jump d[c+1]/d[c].
    dn = dst[:, 1:M]
    gaps = dn[:, 1:] / np.maximum(dn[:, :-1], 1e-6)
    g = np.argmax(gaps, axis=1)  # 0-based gap index
    inner_col = 1 + g  # dst column of the last in-cluster member (d_r)
    outer_col = 2 + g  # dst column of the first non-cluster point (d_{r+1})
    inner = dst[rows, inner_col]
    outer = dst[rows, outer_col]
    jump = outer / np.maximum(inner, 1e-6)
    strong = jump >= JUMP

    # Is the first-excluded neighbour a real same-track cross-image co-observation?
    out_nb = idx[rows, outer_col]
    out_same = (point[out_nb] == point) & pos & (img[out_nb] != img)

    # Same-track cross-image co-observations captured within inner vs outer radius.
    nb_same = (point[idx] == point[:, None]) & pos[:, None] & (img[idx] != img[:, None])
    capt_inner = (nb_same & (dst <= inner[:, None])).sum(1)
    capt_outer = (nb_same & (dst <= outer[:, None])).sum(1)

    return dict(
        name=ws.replace("_ws", ""),
        n=bank.n,
        pos=pos,
        strong=strong,
        d1=dst[:, 1],
        inner=inner,
        outer=outer,
        jump=jump,
        out_same=out_same,
        capt_inner=capt_inner,
        capt_outer=capt_outer,
        T_otsu=float(derive_threshold(dst[:, 1], "otsu")),
    )


def evaluate(d1, pos, T):
    """Separation of in-track vs background if we keep points with d1 <= T."""
    recall = float((d1[pos] <= T).mean())
    bg_drop = float((d1[~pos] > T).mean())
    return recall, bg_drop, (recall + bg_drop) / 2.0


def print_percentiles(r: dict) -> None:
    """Sweep percentiles of the first-excluded d_(r+1) distribution (all points)
    as candidate global radii. Lower percentile -> tighter; for a recall-first
    matcher we care most about keeping recall high (RANSAC filters precision)."""
    d1, pos, outer = r["d1"], r["pos"], r["outer"]
    print(f"  percentile-of-d_(r+1) (all)   global Otsu T={r['T_otsu']:.0f}")
    print(f"  {'pct':>4} {'T':>6} {'recall':>7} {'bg_drop':>8} {'balacc':>7}")
    for p in (10, 20, 30, 40, 50):
        T = float(np.percentile(outer, p))
        rec, drop, bal = evaluate(d1, pos, T)
        print(f"  {p:>3}% {T:>6.1f} {rec:>7.1%} {drop:>8.1%} {bal:>7.3f}")


def print_thresholds(r: dict) -> None:
    """Could the median first-excluded distance serve as the global radius?"""
    d1, pos, s = r["d1"], r["pos"], r["strong"]
    outer = r["outer"]
    cands = [
        ("Otsu(d1)  [current]", r["T_otsu"]),
        ("Youden(d1) [oracle]", youden(d1, pos)[0]),
        ("median d_(r+1) all", float(np.median(outer))),
        ("median d_(r+1) strong", float(np.median(outer[s]))),
        ("median d_(r+1) in-trk", float(np.median(outer[s & pos]))),
    ]
    print(f"  {'threshold':<24} {'T':>6} {'recall':>7} {'bg_drop':>8} {'balacc':>7}")
    for name, T in cands:
        rec, drop, bal = evaluate(d1, pos, T)
        print(f"  {name:<24} {T:>6.1f} {rec:>7.1%} {drop:>8.1%} {bal:>7.3f}")


def print_table(r: dict) -> None:
    sel = r["strong"] & r["pos"]
    inner, outer = r["inner"][sel], r["outer"][sel]
    same = r["out_same"][sel]
    print(f"\n{r['name']}: n={r['n']}  in-track strong cliffs={int(sel.sum())}")
    print(
        f"  median inner d_r={np.median(inner):.0f}  "
        f"first-excluded d_(r+1)={np.median(outer):.0f}  "
        f"jump={np.median(r['jump'][sel]):.2f}  global T={r['T_otsu']:.0f}"
    )
    print(
        f"  first-excluded IS a same-track co-obs: {same.mean():.1%}  "
        f"(=> cliff amputated a real member)"
    )
    print(
        f"  same-track co-obs captured  within inner={r['capt_inner'][sel].mean():.2f}"
        f"  within outer={r['capt_outer'][sel].mean():.2f}"
        f"  (+{r['capt_outer'][sel].mean() - r['capt_inner'][sel].mean():.2f})"
    )


def plot_inner_outer(results, outpath):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, r in zip(axes.ravel(), results):
        sel = r["strong"] & r["pos"]
        inner = r["inner"][sel]
        outer = r["outer"][sel]
        same = r["out_same"][sel]
        hi = np.percentile(outer, 99)
        bins = np.linspace(0, hi, 60)
        ax.hist(
            inner,
            bins=bins,
            density=True,
            alpha=0.55,
            color="C0",
            label="inner radius d_r",
        )
        ax.hist(
            outer[same],
            bins=bins,
            density=True,
            alpha=0.55,
            color="C2",
            label="first-excluded: same-track",
        )
        ax.hist(
            outer[~same],
            bins=bins,
            density=True,
            alpha=0.55,
            color="C1",
            label="first-excluded: background",
        )
        ax.axvline(r["T_otsu"], color="k", ls="--", label=f"global T={r['T_otsu']:.0f}")
        ax.set_title(f"{r['name']}  (cut-too-early {same.mean():.0%})")
        ax.set_xlabel("descriptor distance")
        ax.set_ylabel("density")
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(
        "Inner radius vs first non-cluster point (in-track strong cliffs)", y=0.997
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, dpi=130)
    print(f"\nwrote {outpath}")


def plot_percentiles(results, outpath):
    """recall / background-drop / radius as we sweep the d_(r+1) percentile."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ps = np.arange(5, 91, 5)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, r in zip(axes.ravel(), results):
        d1, pos, outer = r["d1"], r["pos"], r["outer"]
        T = np.percentile(outer, ps)
        rec = np.array([(d1[pos] <= t).mean() for t in T])
        drop = np.array([(d1[~pos] > t).mean() for t in T])
        ax.plot(ps, rec, "o-", color="C0", label="recall (in-track kept)")
        ax.plot(ps, drop, "s-", color="C1", label="background dropped")
        # percentile whose T matches the global Otsu radius
        p_otsu = float(np.interp(r["T_otsu"], T, ps))
        ax.axvline(p_otsu, color="k", ls="--", label=f"Otsu T at p={p_otsu:.0f}")
        ax.axvline(20, color="C3", ls=":", label="p=20")
        ax.set_title(r["name"])
        ax.set_xlabel("percentile of first-excluded d_(r+1)")
        ax.set_ylabel("fraction")
        ax.set_ylim(0, 1.02)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Percentile of first-excluded distance as the global radius", y=0.997)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, dpi=130)
    print(f"\nwrote {outpath}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="accurate")
    ap.add_argument("--datasets", nargs="*", default=DATASETS)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    print(
        "Inner radius d_r vs first non-cluster point d_(r+1), in-track strong "
        "cliffs.\nIf the first-excluded is often a same-track co-obs, the cliff "
        "cuts too early."
    )
    results = [run_one(ws, args.preset) for ws in args.datasets]
    for r in results:
        print_table(r)
        print_thresholds(r)
        print_percentiles(r)
    if len(results) == len(DATASETS):
        out = Path(args.outdir)
        out.mkdir(parents=True, exist_ok=True)
        plot_inner_outer(results, str(out / "exp14_inner_outer.png"))
        plot_percentiles(results, str(out / "exp14_percentiles.png"))


if __name__ == "__main__":
    main()
