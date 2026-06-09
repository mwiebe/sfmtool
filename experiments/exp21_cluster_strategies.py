# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 21 — benchmark local cluster-determination strategies.

Loads the labelled neighbourhood fixtures from exp20 and scores candidate
"given a seed's neighbours, which are cluster members?" rules against the
ground-truth co-observations — fast, no index or solve needed. Add a function to
``STRATEGIES`` and rerun to compare.

A strategy takes the per-seed arrays of one fixture and returns a boolean
``(S, KB)`` predicted-member mask. Scoring (over in-track seeds): precision,
recall, F1 vs ``is_coobs``, plus ``bg_leak`` = average members predicted on
background seeds (should be ~0). Same-image neighbours are never co-observations,
so predicting one is a false positive.

Usage:
    pixi run -e experiments python experiments/exp21_cluster_strategies.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

EPS = 1e-6


def _ranks0(case):
    return np.arange(case["dist"].shape[1])


# --- strategies: case -> (S, KB) bool member mask --------------------------------


def strat_cliff7(case):
    """Largest relative jump among the first 7 neighbours; keep up to it."""
    d = case["dist"][:, :7]
    g = (d[:, 1:] / np.maximum(d[:, :-1], EPS)).argmax(1)  # gap index 0..5
    keep = _ranks0(case)[None, :] <= g[:, None]  # ranks 1..g+1
    return keep & ~case["same_image"]


def strat_radius_cd2(case, c=2.0):
    """Per-point radius = c x the 2nd-nearest distance."""
    d2 = case["dist"][:, 1:2]
    return (case["dist"] <= c * d2) & ~case["same_image"]


def strat_reciprocity(case, r=4):
    """Mutual neighbours: the seed is within the neighbour's r nearest."""
    rr = case["rev_rank"]
    return (rr >= 1) & (rr <= r) & ~case["same_image"]


def strat_recip_and_cd2(case, r=8, c=2.5):
    return strat_reciprocity(case, r) & strat_radius_cd2(case, c)


STRATEGIES = {
    "cliff7": strat_cliff7,
    "c·d2 (c=2)": strat_radius_cd2,
    "recip R≤4": strat_reciprocity,
    "recip&c·d2": strat_recip_and_cd2,
}


def score(case, pred):
    pos = case["is_coobs"]
    in_track = pos.any(1)
    tp = (pred & pos).sum(1).astype(float)
    fp = (pred & ~pos).sum(1).astype(float)
    fn = (~pred & pos).sum(1).astype(float)
    it = in_track
    prec = np.where((tp + fp) > 0, tp / (tp + fp + EPS), 1.0)[it]
    rec = (tp / (tp + fn + EPS))[it]
    f1 = np.where((prec + rec) > 0, 2 * prec * rec / (prec + rec + EPS), 0.0)
    bg_leak = pred[~in_track].sum(1).mean() if (~in_track).any() else 0.0
    return prec.mean(), rec.mean(), f1.mean(), bg_leak


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--casedir", default="cluster_cases")
    args = ap.parse_args()
    files = sorted(Path(args.casedir).glob("*.npz"))
    if not files:
        raise SystemExit(f"no fixtures in {args.casedir}; run exp20 first")
    for fp in files:
        case = dict(np.load(fp, allow_pickle=True))
        S = len(case["seed_global"])
        it = int(case["is_coobs"].any(1).sum())
        print(
            f"\n{case['dataset'].item() if case['dataset'].ndim == 0 else fp.stem}"
            f": {S} seeds ({it} in-track), KB={int(case['kb'])}"
        )
        print(f"  {'strategy':<14} {'prec':>6} {'recall':>7} {'F1':>6} {'bg_leak':>8}")
        for sname, fn in STRATEGIES.items():
            p, r, f, leak = score(case, fn(case))
            print(f"  {sname:<14} {p:>6.3f} {r:>7.3f} {f:>6.3f} {leak:>8.2f}")


if __name__ == "__main__":
    main()
