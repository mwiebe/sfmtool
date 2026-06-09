# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 23 — tune the background-floor strategy on the fixtures.

bgfloor keeps neighbours within `alpha * B`, where `B` is the background scale
estimated as the median distance of the neighbours from rank `b0` onward. This
sweeps `alpha` and `b0` over the exp20 fixtures (instant — no index/solve),
ranks configs by mean F1 across datasets, and shows the per-dataset breakdown for
the default, the best config, and a reciprocity-gated variant.

Usage:
    pixi run -e experiments python experiments/exp23_bgfloor_sweep.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from exp21_cluster_strategies import _cross, score

ALPHAS = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
B0S = [8, 16, 24, 32]


def bgfloor(case, alpha, b0):
    b = np.median(case["dist"][:, b0:], axis=1, keepdims=True)
    return (case["dist"] <= alpha * b) & _cross(case)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--casedir", default="cluster_cases")
    args = ap.parse_args()
    cases = {
        fp.stem: dict(np.load(fp, allow_pickle=True))
        for fp in sorted(Path(args.casedir).glob("*.npz"))
    }
    names = list(cases)

    def mean_f1(alpha, b0):
        return float(
            np.mean([score(c, bgfloor(c, alpha, b0))[2] for c in cases.values()])
        )

    grid = sorted(((mean_f1(a, b0), a, b0) for a in ALPHAS for b0 in B0S), reverse=True)
    print("top configs by mean F1 across datasets (alpha, b0):")
    for f, a, b0 in grid[:8]:
        print(f"  alpha={a:.1f} b0={b0:<2}  meanF1={f:.3f}")

    best_f, best_a, best_b0 = grid[0]

    def breakdown(label, fn):
        print(f"\n{label}:")
        print(f"  {'dataset':<18} {'prec':>6} {'recall':>7} {'F1':>6} {'bg_leak':>8}")
        for n in names:
            p, r, f, leak = score(cases[n], fn(cases[n]))
            print(f"  {n:<18} {p:>6.3f} {r:>7.3f} {f:>6.3f} {leak:>8.2f}")

    breakdown("default (alpha=0.7, b0=24)", lambda c: bgfloor(c, 0.7, 24))
    breakdown(
        f"best (alpha={best_a:.1f}, b0={best_b0})",
        lambda c: bgfloor(c, best_a, best_b0),
    )
    breakdown(
        "best & reciprocity (rev_rank≥1)",
        lambda c: bgfloor(c, best_a, best_b0) & (c["rev_rank"] >= 1),
    )


if __name__ == "__main__":
    main()
