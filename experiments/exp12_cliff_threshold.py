# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 12 — derive the match radius from each point's neighbour "cliff".

Instead of fitting Otsu/GMM on the global d1 distribution, look at every point's
own sorted neighbours and find the *cliff*: the rank where the distance suddenly
jumps from "near" (its true co-observations) to "far" (background). The distance
just before that jump is a per-point estimate of the local cluster radius. We
aggregate those per-point radii (P50/P75/P90) into one global T and compare it,
and the in-track/background separation it gives, to the current Otsu(d1) radius
and the Youden-optimal d1 oracle.

Cost note: the matcher already queries k=17 (exp03), so reading 8 neighbours and
scanning 7 gaps is a free by-product — no extra index work, no EM fit.

Usage:
    pixi run -e experiments python experiments/exp12_cliff_threshold.py
"""

from __future__ import annotations

import argparse
import glob

import numpy as np

from exp05_cluster_match import derive_threshold
from exp11_ratio_sweep import youden
from sfm_descriptors import load_descriptor_bank

DATASETS = [
    "seoul_bull_ws",
    "seattle_backyard_ws",
    "kerry_park_ws",
    "dino_dog_toy_ws",
]


def cliff(dn: np.ndarray):
    """Per-point cliff in sorted neighbour distances ``dn`` (N, m), ascending.

    Returns ``(cliff_dist, cliff_rank, jump)``: the distance just before the
    largest *relative* jump d[r+1]/d[r], that rank r (1-based count of near
    members), and the jump ratio there (cliff strength).
    """
    denom = np.maximum(dn[:, :-1], 1e-6)
    gaps = dn[:, 1:] / denom  # (N, m-1) ratio between consecutive neighbours
    r = np.argmax(gaps, axis=1)  # 0-based index of the gap with the biggest jump
    rows = np.arange(len(dn))
    cliff_dist = dn[rows, r]
    jump = gaps[rows, r]
    return cliff_dist, r + 1, jump


def evaluate(d1: np.ndarray, pos: np.ndarray, T: float):
    """Separation of in-track vs background if we keep points with d1 <= T."""
    recall = float((d1[pos] <= T).mean())
    bg_drop = float((d1[~pos] > T).mean())
    return recall, bg_drop, (recall + bg_drop) / 2.0


def run_one(ws: str, preset: str) -> None:
    from sfmtool import KdForest

    path = sorted(glob.glob(f"../{ws}/sfmr/*solve*.sfmr"))[0]
    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    _, dst = KdForest(desc, preset=preset).query(desc, k=8)
    pos = bank.point_label >= 0
    d1 = dst[:, 1]
    dn = dst[:, 1:8]  # 7 nearest *others*, ascending

    cdist, crank, jump = cliff(dn)
    # Unsupervised gate: a point "has a cluster" only if its cliff is a real
    # drop (>=30% jump). Background points have near-flat neighbours (jump ~1).
    strong = jump >= 1.3

    # Candidate global radii.
    cands = [
        ("Otsu(d1)  [current]", derive_threshold(d1, "otsu")),
        ("Youden(d1) [oracle]", youden(d1, pos)[0]),
        ("cliff P90 (all)", float(np.percentile(cdist, 90))),
        ("cliff P90 (strong)", float(np.percentile(cdist[strong], 90))),
        ("cliff P75 (strong)", float(np.percentile(cdist[strong], 75))),
        ("cliff P90 (in-trk)", float(np.percentile(cdist[pos], 90))),
    ]

    print(
        f"\n{ws.replace('_ws', '')}: n={bank.n}  in-track={pos.mean():.1%}  "
        f"strong-cliff={strong.mean():.1%} (of which in-track={pos[strong].mean():.1%})"
    )
    print(
        f"  cliff rank  median in-track={np.median(crank[pos]):.0f}  "
        f"background={np.median(crank[~pos]):.0f}   "
        f"cliff dist  median in-track={np.median(cdist[pos]):.0f}  "
        f"background={np.median(cdist[~pos]):.0f}"
    )
    print(f"  {'method':<22} {'T':>6} {'recall':>7} {'bg_drop':>8} {'balacc':>7}")
    for name, T in cands:
        rec, drop, bal = evaluate(d1, pos, T)
        print(f"  {name:<22} {T:>6.1f} {rec:>7.1%} {drop:>8.1%} {bal:>7.3f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="accurate")
    ap.add_argument("--datasets", nargs="*", default=DATASETS)
    args = ap.parse_args()
    print(
        "Radius from per-point neighbour cliff vs the current Otsu(d1) radius.\n"
        "recall = in-track kept (d1<=T); bg_drop = background rejected;\n"
        "balacc = (recall+bg_drop)/2 at that single radius."
    )
    for ws in args.datasets:
        run_one(ws, args.preset)


if __name__ == "__main__":
    main()
