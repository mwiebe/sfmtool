# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 11 — nearest-vs-Nth ratio sweep (N = 2..7).

The POC's core matcher does *not* use a ratio at all: it thresholds the absolute
nearest-other distance `d1` at a data-derived radius `T`. The only ratio in the
codebase is the optional `--prefilter` (`d1/d5 > 0.85` drops isolated points) —
nearest-vs-5th, not the classic Lowe nearest-vs-2nd.

This sweep asks how the discriminative power of the ratio `r_N = d1 / d_N` shifts
as N goes 2 -> 7, using the solve's tracks as ground truth (in-track vs the
background the solve discarded). For each N we report how separable in-track and
background are under r_N (ROC-AUC), the median ratio for each class, and the
Youden-optimal ratio threshold with its recall / background-rejection. The
absolute-`d1` radius (what the matcher actually uses) is the baseline row.

Usage:
    pixi run -e experiments python experiments/exp11_ratio_sweep.py
"""

from __future__ import annotations

import argparse
import glob
from time import perf_counter

import numpy as np

from exp05_cluster_match import K
from sfm_descriptors import load_descriptor_bank

DATASETS = [
    "seoul_bull_ws",
    "seattle_backyard_ws",
    "kerry_park_ws",
    "dino_dog_toy_ws",
]


def auc(score: np.ndarray, pos: np.ndarray) -> float:
    """ROC-AUC = P(score[pos] > score[neg]) via average ranks (Mann-Whitney)."""
    order = np.argsort(score, kind="stable")
    ranks = np.empty(len(score), dtype=np.float64)
    ranks[order] = np.arange(1, len(score) + 1)
    # average ranks within tied score groups
    s = score[order]
    ties = np.r_[0, np.flatnonzero(np.diff(s)) + 1, len(s)]
    for a, b in zip(ties[:-1], ties[1:]):
        if b - a > 1:
            ranks[order[a:b]] = (a + 1 + b) / 2.0
    npos = int(pos.sum())
    nneg = len(score) - npos
    if npos == 0 or nneg == 0:
        return float("nan")
    return (ranks[pos].sum() - npos * (npos + 1) / 2.0) / (npos * nneg)


def youden(value: np.ndarray, pos: np.ndarray):
    """Best 'in-track iff value <= thr' split: returns (thr, recall, bg_drop)."""
    order = np.argsort(value, kind="stable")
    v = value[order]
    p = pos[order].astype(np.int64)
    tp = np.cumsum(p)
    fp = np.cumsum(1 - p)
    P, Nn = tp[-1], fp[-1]
    tpr = tp / max(P, 1)
    fpr = fp / max(Nn, 1)
    j = tpr - fpr
    i = int(np.argmax(j))
    return float(v[i]), float(tpr[i]), float(1.0 - fpr[i])


def time_query(forest, q: np.ndarray, k: int, repeats: int) -> float:
    """Median wall-clock seconds for one full ``query(q, k)`` over the forest."""
    ts = []
    for _ in range(repeats):
        t0 = perf_counter()
        forest.query(q, k=k)
        ts.append(perf_counter() - t0)
    return float(np.median(ts))


def run_one(ws: str, preset: str, repeats: int) -> dict:
    """Sweep d1/dN separability and per-k query cost for one dataset."""
    from sfmtool import KdForest

    path = sorted(glob.glob(f"../{ws}/sfmr/*solve*.sfmr"))[0]
    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    forest = KdForest(desc, preset=preset)
    # k must reach the 7th nearest *other* (col 0 is self).
    idx, dst = forest.query(desc, k=max(K, 8))
    pos = bank.point_label >= 0
    d1 = dst[:, 1]

    # Cost along the sweep axis: to use d1/dN you must retrieve k = N+1
    # neighbours; the absolute-d1 matcher only needs k = 2.
    ns = list(range(2, 8))
    ks = [2, *(n + 1 for n in ns)]  # k=2 for abs-d1; k=N+1 to reach dN
    t_us = {k: time_query(forest, desc, k, repeats) / bank.n * 1e6 for k in ks}

    rows = []  # (label, k, auc, intrk_md, bg_md, thr, recall, bg_drop, us_per_q)
    a = auc(-d1, pos)
    thr, rec, drop = youden(d1, pos)
    rows.append(
        (
            "d1 (abs)",
            2,
            a,
            np.median(d1[pos]),
            np.median(d1[~pos]),
            thr,
            rec,
            drop,
            t_us[2],
        )
    )
    for n in ns:
        r = d1 / np.maximum(dst[:, n], 1e-6)
        a = auc(-r, pos)
        thr, rec, drop = youden(r, pos)
        rows.append(
            (
                f"d1/d{n}",
                n + 1,
                a,
                np.median(r[pos]),
                np.median(r[~pos]),
                thr,
                rec,
                drop,
                t_us[n + 1],
            )
        )

    name = ws.replace("_ws", "")
    print(f"\n{name}: n={bank.n}  in-track={pos.mean():.1%}")
    print(
        f"  {'signal':<12} {'k':>2} {'AUC':>6} {'intrk_md':>9} {'bg_md':>7} "
        f"{'thr':>6} {'recall':>7} {'bg_drop':>8} {'us/query':>9}"
    )
    for lbl, k, a, im, bm, thr, rec, drop, us in rows:
        fmt = ">9.1f" if lbl == "d1 (abs)" else ">9.3f"
        bfmt = ">7.1f" if lbl == "d1 (abs)" else ">7.3f"
        tfmt = ">6.1f" if lbl == "d1 (abs)" else ">6.3f"
        print(
            f"  {lbl:<12} {k:>2} {a:>6.3f} {im:{fmt}} {bm:{bfmt}} "
            f"{thr:{tfmt}} {rec:>7.1%} {drop:>8.1%} {us:>9.2f}"
        )

    return dict(
        name=name,
        ns=ns,
        auc_abs=rows[0][2],
        us_abs=t_us[2],
        auc=[r[2] for r in rows[1:]],
        us=[t_us[n + 1] for n in ns],
    )


def plot_sweep(results: list[dict], outpath: str) -> None:
    """Per-dataset AUC (separability) and per-query cost vs N, on twin axes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    handles = None
    for i, (ax, r) in enumerate(zip(axes.ravel(), results)):
        ax.plot(r["ns"], r["auc"], "o-", color="C0", label="d1/dN AUC")
        ax.axhline(r["auc_abs"], ls="--", color="C2", label="d1 abs AUC (k=2)")
        ax.set_ylabel("AUC (separability)", color="C0")
        ax.tick_params(axis="y", labelcolor="C0")
        ax.set_xlabel("N  (ratio d1/dN; query k = N+1)")
        ax.set_title(r["name"])
        ax.grid(alpha=0.3)
        axc = ax.twinx()
        axc.plot(r["ns"], r["us"], "s-", color="C3", label="query cost (µs)")
        axc.axhline(r["us_abs"], ls=":", color="C3", alpha=0.6)
        axc.set_ylabel("µs / query", color="C3")
        axc.tick_params(axis="y", labelcolor="C3")
        axc.set_ylim(bottom=0)
        if i == 0:
            h0, l0 = ax.get_legend_handles_labels()
            h1, l1 = axc.get_legend_handles_labels()
            handles = (h0 + h1, l0 + l1)
    fig.legend(*handles, loc="upper center", ncol=3, frameon=False)
    fig.suptitle("Nearest-vs-Nth ratio: separability vs query cost", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(outpath, dpi=130)
    print(f"\nwrote {outpath}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="accurate")
    ap.add_argument("--datasets", nargs="*", default=DATASETS)
    ap.add_argument(
        "--repeats", type=int, default=5, help="query-timing repeats (median reported)"
    )
    ap.add_argument(
        "--plot",
        default="out/exp11_ratio_sweep.png",
        help="output figure path (empty string to skip)",
    )
    args = ap.parse_args()
    print(
        "Separability of in-track vs background descriptors, and query cost.\n"
        "AUC = P(in-track scored above background); higher = more separable.\n"
        "thr = Youden-optimal 'in-track iff signal <= thr'; recall = in-track\n"
        "kept; bg_drop = background rejected; us/query = index query at k=N+1."
    )
    results = [run_one(ws, args.preset, args.repeats) for ws in args.datasets]
    if args.plot and len(results) == len(DATASETS):
        from pathlib import Path

        Path(args.plot).parent.mkdir(parents=True, exist_ok=True)
        plot_sweep(results, args.plot)


if __name__ == "__main__":
    main()
