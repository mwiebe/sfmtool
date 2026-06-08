# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 07 — how do materialized clusters line up with the solve's tracks?

Builds the in-tree KdForest index, derives the radius `T` from the data (Otsu),
materialises clusters with density-ordered seeding (no transitive merge),
optionally mean-shift refined, and scores those clusters against the solve's
ground-truth tracks — purity,
false-merge rate, one-feature-per-image rate, track recovery, fragmentation, and
coverage — for all four datasets.

Usage:
    pixi run -e experiments python experiments/exp07_cluster_vs_tracks.py
"""

from __future__ import annotations

import argparse
import glob

import numpy as np

from exp03_radius_clusters import score
from exp05_cluster_match import K, cliff_threshold, derive_threshold
from seed_cluster import seed_claim_clusters
from sfm_descriptors import load_descriptor_bank

DATASETS = [
    "seoul_bull_ws",
    "seattle_backyard_ws",
    "kerry_park_ws",
    "dino_dog_toy_ws",
]


def run_one(ws: str, threshold: str, t_scale: float, preset: str, refine: int):
    from sfmtool import KdForest

    path = sorted(glob.glob(f"../{ws}/sfmr/*solve*.sfmr"))[0]
    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    forest = KdForest(desc, preset=preset)
    idx, dst = forest.query(desc, k=K)
    idx = idx.astype(np.int64)

    if threshold == "cliff":
        T = cliff_threshold(dst) * t_scale
    else:
        T = derive_threshold(dst[:, 1], threshold) * t_scale
    labels = seed_claim_clusters(
        idx, dst, bank.image_label, T, descriptors=desc,
        refine_iters=refine, forest=forest,
    )

    # Reuse exp03's scorer: give every unclustered descriptor its own singleton
    # label so cluster sizes are correct, then score the multi-member clusters.
    n = bank.n
    comp = labels.copy()
    unl = np.flatnonzero(labels < 0)
    next_id = int(labels.max()) + 1 if (labels >= 0).any() else 0
    comp[unl] = np.arange(len(unl), dtype=np.int64) + next_id
    sizes = np.bincount(comp, minlength=n)
    res = score(bank, comp, sizes, np.ones(n, dtype=bool), T)

    n_tracks = len(np.unique(bank.point_label[bank.point_label >= 0]))
    res.update(name=ws.replace("_ws", ""), n=n, T=T, tracks=n_tracks)
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--threshold", default="cliff",
                    help="radius: cliff (p50; default) / otsu / gmm / float")
    ap.add_argument("--t-scale", type=float, default=1.0)
    ap.add_argument("--preset", default="accurate")
    ap.add_argument("--refine", type=int, default=0,
                    help="mean-shift re-query iterations (0 = seed only)")
    ap.add_argument("--datasets", nargs="*", default=DATASETS)
    args = ap.parse_args()

    rows = []
    for ws in args.datasets:
        print(f"clustering {ws} (refine={args.refine}) ...", flush=True)
        rows.append(
            run_one(ws, args.threshold, args.t_scale, args.preset, args.refine)
        )

    print(f"\nClusters vs solve tracks ({args.threshold} T × {args.t_scale}, "
          f"refine={args.refine}):\n")
    hdr = (f"  {'dataset':<16} {'T':>5} {'clusters':>8} {'tracks':>6} "
           f"{'purity':>7} {'falseM':>7} {'1/img':>6} {'recov':>6} "
           f"{'frags':>6} {'cover':>6}")
    print(hdr)
    for r in rows:
        print(f"  {r['name']:<16} {r['T']:>5.0f} {r['components']:>8} "
              f"{r['tracks']:>6} {r['dom_cov']:>7.3f} {r['false_merge_frac']:>7.1%} "
              f"{r['img_unique_frac']:>6.1%} {r['recovery']:>6.3f} "
              f"{r['fragments_per_track']:>6.2f} {r['in_track_coverage']:>6.1%}")
    print("\npurity = in-track members in their cluster's dominant track; "
          "1/img = clusters with ≤1 feature per image;\nrecov = mean track "
          "completeness; frags = clusters a track spans; cover = in-track "
          "descriptors landed in a cluster.")


if __name__ == "__main__":
    main()
