# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 07 — how do materialized clusters line up with the solve's tracks?

Builds the in-tree KdForest index, materialises clusters by density-ordered
seeding (no transitive merge), and scores those clusters against the solve's
ground-truth tracks for all four datasets. Two membership rules:

  - ``bgfloor`` (default): each seed keeps neighbours within its own background
    floor ``alpha * dist[d]`` (d=28, alpha=0.8) — the production rule.
  - ``global``: a single data-derived radius (cliff p50) for the whole corpus.

Reports, per dataset: the number of solve tracks, the number of clusters we form,
how many clusters span >1 track (overlap), how many contain no track member
(non-track), and mean track recovery.

Usage:
    pixi run -e experiments python experiments/exp07_cluster_vs_tracks.py
"""

from __future__ import annotations

import argparse
import glob

import numpy as np

from exp03_radius_clusters import score
from exp05_cluster_match import cliff_threshold
from seed_cluster import seed_claim_clusters
from sfm_descriptors import load_descriptor_bank
from sfmtool._sfmtool import background_floor_clusters

DATASETS = [
    "seoul_bull_ws",
    "seattle_backyard_ws",
    "kerry_park_ws",
    "dino_dog_toy_ws",
]

D_RANK = 28  # background rank: B_i = dist[i, 28] (production default in main is 10)
BG_ALPHA = 0.8


def production_labels(bank, preset: str) -> np.ndarray:
    """Per-descriptor cluster labels from the production background-floor matcher.

    Runs the shipped ``sfmtool.background_floor_clusters`` (Rust) and scatters its
    CSR clusters back onto the bank's concatenated rows: ``labels[row] = cluster``
    for every member, ``-1`` for descriptors no cluster claimed.
    """
    corpus = np.ascontiguousarray(bank.descriptors)
    image_starts = bank.image_starts
    cluster_starts, member_images, member_features = background_floor_clusters(
        corpus, image_starts, d=D_RANK, alpha=BG_ALPHA, min_size=2, preset=preset
    )
    rows = image_starts[member_images].astype(np.int64) + member_features
    counts = np.diff(cluster_starts.astype(np.int64))
    member_cluster = np.repeat(np.arange(len(counts), dtype=np.int64), counts)
    labels = np.full(bank.n, -1, dtype=np.int64)
    labels[rows] = member_cluster
    return labels


def run_one(ws: str, mode: str, preset: str):
    from sfmtool import KdForest

    path = sorted(glob.glob(f"../{ws}/sfmr/*solve*.sfmr"))[0]
    bank = load_descriptor_bank(path)

    if mode == "bgfloor":
        # The production rule: cluster with the shipped Rust matcher itself.
        labels = production_labels(bank, preset)
    else:
        # The global-threshold alternative we considered (POC numpy clustering).
        desc = np.ascontiguousarray(bank.descriptors)
        forest = KdForest(desc, preset=preset)
        idx, dst = forest.query(desc, k=17)
        idx = idx.astype(np.int64)
        labels = seed_claim_clusters(idx, dst, bank.image_label, cliff_threshold(dst))

    # Singletons for unclustered descriptors so cluster sizes are correct.
    comp = labels.copy()
    unl = np.flatnonzero(labels < 0)
    next_id = int(labels.max()) + 1 if (labels >= 0).any() else 0
    comp[unl] = np.arange(len(unl), dtype=np.int64) + next_id
    sizes = np.bincount(comp, minlength=bank.n)
    res = score(bank, comp, sizes, np.ones(bank.n, dtype=bool), 0.0)

    n_clusters = int(labels.max()) + 1 if (labels >= 0).any() else 0
    n_tracks = len(np.unique(bank.point_label[bank.point_label >= 0]))
    return {
        "name": ws.replace("_ws", ""),
        "tracks": n_tracks,
        "clusters": n_clusters,
        "overlaps": res.get("false_merge_comps", 0),
        "non_track": n_clusters - res.get("components", 0),
        "recovery": res.get("recovery", 0.0),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=["bgfloor", "global"], default="bgfloor")
    ap.add_argument("--preset", default="accurate")
    ap.add_argument("--datasets", nargs="*", default=DATASETS)
    args = ap.parse_args()

    rows = []
    for ws in args.datasets:
        print(f"clustering {ws} (mode={args.mode}) ...", flush=True)
        rows.append(run_one(ws, args.mode, args.preset))

    print(f"\nClusters vs solve tracks (mode={args.mode}):\n")
    print(
        f"  {'dataset':<16} {'tracks':>7} {'clusters':>8} {'overlaps':>8} "
        f"{'non-track':>9} {'recovery':>8}"
    )
    for r in rows:
        print(
            f"  {r['name']:<16} {r['tracks']:>7} {r['clusters']:>8} "
            f"{r['overlaps']:>8} {r['non_track']:>9} {r['recovery']:>8.2f}"
        )


if __name__ == "__main__":
    main()
