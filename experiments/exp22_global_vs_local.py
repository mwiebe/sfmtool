# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 22 — global threshold vs local determination, for the leaders.

The two leading local strategies (background-floor and absolute-gap cliff,
exp21) decide membership *per seed*. This asks two things on the labelled
fixtures (exp20), scored against each solve's tracks:

  1. How to turn a local rule into a single global radius `T`: take each seed's
     effective radius (the farthest member it keeps) and aggregate (median) over
     in-track seeds. We then score that global `T` and compare it to the local
     rule it came from, and to the *best possible* global `T` (grid-searched to
     maximise F1 — the ceiling of any single threshold).
  2. Whether local beats global at all: compare the local rules' F1 to the best
     global F1.

Usage:
    pixi run -e experiments python experiments/exp22_global_vs_local.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from exp21_cluster_strategies import (
    _cross,
    score,
    strat_bgfloor,
    strat_cliff_abs,
)

LEADERS = {"bgfloor": strat_bgfloor, "cliff_abs": strat_cliff_abs}


def per_seed_radius(case, pred):
    """Effective radius each seed used = farthest member it kept (nan if none)."""
    d = np.where(pred, case["dist"], -np.inf)
    r = d.max(1)
    r[~np.isfinite(r)] = np.nan
    return r


def global_from_local(case, pred):
    """Median effective radius over in-track seeds → one global T."""
    it = case["is_coobs"].any(1)
    r = per_seed_radius(case, pred)[it]
    r = r[np.isfinite(r)]
    return float(np.median(r)) if r.size else 0.0


def apply_global(case, t):
    return (case["dist"] <= t) & _cross(case)


def best_global(case):
    """Grid-search the single global T that maximises mean F1."""
    d = case["dist"][_cross(case)]
    cand = np.quantile(d, np.linspace(0.02, 0.7, 80))
    best_f, best_t = -1.0, 0.0
    for t in cand:
        _, _, f, _ = score(case, apply_global(case, t))
        if f > best_f:
            best_f, best_t = f, float(t)
    return best_t, best_f


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--casedir", default="cluster_cases")
    args = ap.parse_args()
    for fp in sorted(Path(args.casedir).glob("*.npz")):
        case = dict(np.load(fp, allow_pickle=True))
        print(f"\n{fp.stem}:")
        print(f"  {'method':<22} {'T':>6} {'prec':>6} {'recall':>7} {'F1':>6}")
        for lname, fn in LEADERS.items():
            pred = fn(case)
            p, r, f, _ = score(case, pred)
            print(f"  {lname + ' (local)':<22} {'—':>6} {p:>6.3f} {r:>7.3f} {f:>6.3f}")
            t = global_from_local(case, pred)
            pg, rg, fg, _ = score(case, apply_global(case, t))
            print(
                f"  {'  → global median T':<22} {t:>6.0f} {pg:>6.3f} {rg:>7.3f} "
                f"{fg:>6.3f}"
            )
        bt, _ = best_global(case)
        pb, rb, fb, _ = score(case, apply_global(case, bt))
        print(
            f"  {'best global T (oracle)':<22} {bt:>6.0f} {pb:>6.3f} {rb:>7.3f} "
            f"{fb:>6.3f}"
        )


if __name__ == "__main__":
    main()
