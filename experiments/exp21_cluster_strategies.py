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
from functools import partial
from pathlib import Path

import numpy as np

EPS = 1e-6


def _cross(case):
    """Cross-image neighbours — the only ones that can be members."""
    return ~case["same_image"]


def _keep_upto(case, cut):
    """Members = neighbours at rank ≤ cut (0-based, inclusive)."""
    r = np.arange(case["dist"].shape[1])[None, :]
    return r <= cut[:, None]


# --- strategies: case -> (S, KB) bool member mask --------------------------------
# Each may use only dist / image / same_image / rev_rank (never point/is_coobs).


# -- gap / cliff on the sorted distance profile --
def strat_cliff(case, w=7):
    """Cut at the largest *relative* jump d[i+1]/d[i] within the first w."""
    d = case["dist"][:, :w]
    g = (d[:, 1:] / np.maximum(d[:, :-1], EPS)).argmax(1)
    return _keep_upto(case, g) & _cross(case)


def strat_cliff_abs(case, w=16):
    """Cut at the largest *absolute* gap d[i+1]-d[i] within the first w."""
    d = case["dist"][:, :w]
    g = (d[:, 1:] - d[:, :-1]).argmax(1)
    return _keep_upto(case, g) & _cross(case)


def strat_outermost(case, w=16, tau=1.3):
    """Cut at the *last* jump (within w) whose ratio ≥ tau — skip intra-track
    gaps to the final drop into background; if none, keep only the nearest."""
    d = case["dist"][:, :w]
    ratios = d[:, 1:] / np.maximum(d[:, :-1], EPS)
    cols = np.arange(ratios.shape[1])[None, :]
    g = np.where(ratios >= tau, cols, -1).max(1)
    return _keep_upto(case, np.maximum(g, 0)) & _cross(case)


def strat_knee(case, w=16):
    """Kneedle: cut at the rank of max gap between the distance curve and the
    chord from the 1st to the w-th neighbour (the bend into the plateau)."""
    d = case["dist"][:, :w]
    x = np.arange(w) / (w - 1)
    chord = d[:, :1] + (d[:, -1:] - d[:, :1]) * x[None, :]
    g = (chord - d).argmax(1)
    return _keep_upto(case, g) & _cross(case)


# -- radius relative to a per-point scale --
def strat_radius_cd2(case, c=2.0):
    """Keep neighbours within c x the 2nd-nearest distance."""
    return (case["dist"] <= c * case["dist"][:, 1:2]) & _cross(case)


def strat_radius_d1(case, tau=1.5):
    """Keep neighbours within tau x the nearest distance (Lowe-flavoured)."""
    return (case["dist"] <= tau * case["dist"][:, :1]) & _cross(case)


def strat_bgfloor(case, alpha=0.8, b0=8):
    """Estimate background scale B = median distance of the neighbours from rank
    b0 onward; keep neighbours within alpha·B. Tuned defaults (exp23)."""
    b = np.median(case["dist"][:, b0:], axis=1, keepdims=True)
    return (case["dist"] <= alpha * b) & _cross(case)


def strat_topk(case, k=4):
    """Keep the k nearest cross-image neighbours (fixed multiplicity)."""
    cr = _cross(case)
    return cr & (np.cumsum(cr, axis=1) <= k)


# -- reciprocity / mutual nearest-neighbour --
def strat_recip(case, r=4):
    """Mutual: the seed is within the neighbour's r nearest."""
    rr = case["rev_rank"]
    return (rr >= 1) & (rr <= r) & _cross(case)


def strat_mutual_best(case):
    """Strict mutual nearest: each is among the other's very nearest."""
    return (case["rev_rank"] == 1) & _cross(case)


# -- combinations (intersect a recall-y radius with a precise gate) --
def strat_recip_cd2(case, r=8, c=2.5):
    return strat_recip(case, r) & strat_radius_cd2(case, c)


def strat_outer_recip(case):
    return strat_outermost(case) & (case["rev_rank"] >= 1)


def strat_knee_recip(case):
    return strat_knee(case) & (case["rev_rank"] >= 1)


def strat_cliffabs_recip(case):
    return strat_cliff_abs(case) & (case["rev_rank"] >= 1)


STRATEGIES = {
    "cliff7": strat_cliff,
    "cliff16": partial(strat_cliff, w=16),
    "cliff_abs16": strat_cliff_abs,
    "outermost": strat_outermost,
    "knee": strat_knee,
    "c·d2 (2.0)": strat_radius_cd2,
    "c·d2 (1.5)": partial(strat_radius_cd2, c=1.5),
    "tau·d1 (1.5)": strat_radius_d1,
    "bgfloor": strat_bgfloor,
    "top4": strat_topk,
    "mutual-best": strat_mutual_best,
    "recip R≤4": strat_recip,
    "recip R≤8": partial(strat_recip, r=8),
    "recip&c·d2": strat_recip_cd2,
    "outer&recip": strat_outer_recip,
    "knee&recip": strat_knee_recip,
    "cliffabs&recip": strat_cliffabs_recip,
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
