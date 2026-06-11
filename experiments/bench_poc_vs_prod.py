# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Benchmark the POC numpy matcher vs the merged Rust implementation.

Both compute the same thing — background-floor track clusters expanded to
per-image-pair matches — over one descriptor corpus, at the same parameters:

  POC   sfmtool.KdForest index (Rust) + numpy clustering (seed_claim_clusters)
        + numpy C(m,2) pair expansion.
  PROD  sfmtool.background_floor_clusters + clusters_to_pair_matches
        (all in crates/sfmtool-core/src/cluster_match, index build included).

Reports per dataset: wall-clock time (broken into index / cluster / pairs where
each side allows it) and the output sizes (clusters, image pairs, matches) so the
two implementations can be checked for behavioural equivalence.

Usage:
    pixi run -e experiments python experiments/bench_poc_vs_prod.py
    pixi run -e experiments python experiments/bench_poc_vs_prod.py \
        --datasets seoul_bull_ws --d 28 --alpha 0.8
"""

from __future__ import annotations

import argparse
import glob
from time import perf_counter

import numpy as np

from seed_cluster import seed_claim_clusters
from sfm_descriptors import load_descriptor_bank
from sfmtool._sfmtool import background_floor_clusters, clusters_to_pair_matches

DATASETS = [
    "seoul_bull_ws",
    "seattle_backyard_ws",
    "kerry_park_ws",
    "dino_dog_toy_ws",
]


def expand_pairs(bank, labels: np.ndarray, min_size: int):
    """numpy C(m,2) cross-image pair expansion of a label vector — the inner loop
    of exp05.build_cluster_matches_arrays, isolated so it can be timed on its own.
    Returns (n_clusters, n_pairs, n_matches)."""
    img = bank.image_label
    valid = np.flatnonzero(labels >= 0)
    if valid.size == 0:
        return 0, 0, 0
    order = valid[np.argsort(labels[valid], kind="stable")]
    clab = labels[order]
    bounds = np.r_[0, np.flatnonzero(np.diff(clab)) + 1, len(order)]
    pair_set: dict[tuple[int, int], int] = {}
    n_clusters = 0
    for s, e in zip(bounds[:-1], bounds[1:]):
        members = order[s:e]
        if len(members) < min_size:
            continue
        n_clusters += 1
        imgs = img[members]
        o = np.argsort(imgs, kind="stable")
        imgs = imgs[o]
        for a in range(len(imgs)):
            for b in range(a + 1, len(imgs)):
                ia, ib = int(imgs[a]), int(imgs[b])
                if ia != ib:
                    pair_set[(ia, ib)] = pair_set.get((ia, ib), 0) + 1
    return n_clusters, len(pair_set), sum(pair_set.values())


def run_poc(bank, d: int, alpha: float, preset: str, min_size: int):
    from sfmtool import KdForest

    desc = np.ascontiguousarray(bank.descriptors)
    t0 = perf_counter()
    forest = KdForest(desc, preset=preset)
    idx, dst = forest.query(desc, k=d + 1)
    idx = idx.astype(np.int64)
    t_index = perf_counter() - t0

    radius = (alpha * dst[:, d])[:, None]  # per-point background-floor radius
    t1 = perf_counter()
    labels = seed_claim_clusters(idx, dst, bank.image_label, radius, min_size=min_size)
    t_cluster = perf_counter() - t1

    t2 = perf_counter()
    n_clusters, n_pairs, n_matches = expand_pairs(bank, labels, min_size)
    t_pairs = perf_counter() - t2

    return dict(
        index=t_index,
        cluster=t_cluster,
        pairs=t_pairs,
        total=t_index + t_cluster + t_pairs,
        clusters=n_clusters,
        n_pairs=n_pairs,
        matches=n_matches,
    )


def run_prod(bank, d: int, alpha: float, preset: str, min_size: int):
    corpus = np.ascontiguousarray(bank.descriptors)
    image_starts = bank.image_starts

    t0 = perf_counter()
    cluster_starts, member_images, member_features = background_floor_clusters(
        corpus, image_starts, d=d, alpha=alpha, min_size=min_size, preset=preset
    )
    t_cluster = perf_counter() - t0  # index build + k-NN + clustering, all in Rust

    t1 = perf_counter()
    ii, mc, _, _ = clusters_to_pair_matches(
        cluster_starts, member_images, member_features, corpus, image_starts
    )
    t_pairs = perf_counter() - t1

    return dict(
        cluster=t_cluster,
        pairs=t_pairs,
        total=t_cluster + t_pairs,
        clusters=len(cluster_starts) - 1,
        n_pairs=len(ii),
        matches=int(mc.sum()),
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--datasets", nargs="*", default=DATASETS)
    ap.add_argument("--d", type=int, default=28)
    ap.add_argument("--alpha", type=float, default=0.8)
    ap.add_argument("--min-size", type=int, default=2)
    ap.add_argument("--preset", default="accurate")
    args = ap.parse_args()

    rows = []
    for ws in args.datasets:
        path = sorted(glob.glob(f"../{ws}/sfmr/*solve*.sfmr"))[0]
        bank = load_descriptor_bank(path)
        print(f"benchmarking {ws} (n={bank.n}) ...", flush=True)
        poc = run_poc(bank, args.d, args.alpha, args.preset, args.min_size)
        prod = run_prod(bank, args.d, args.alpha, args.preset, args.min_size)
        rows.append((ws.replace("_ws", ""), bank.n, poc, prod))

    print(f"\nPOC (numpy) vs PROD (Rust) — d={args.d} alpha={args.alpha}\n")
    print(
        f"  {'dataset':<14} {'descriptors':>11} | "
        f"{'POC total':>9} {'PROD total':>10} {'speedup':>7} | "
        f"{'clusters P/R':>14} {'matches P/R':>16}"
    )
    for name, n, poc, prod in rows:
        speedup = poc["total"] / prod["total"] if prod["total"] else float("nan")
        print(
            f"  {name:<14} {n:>11,} | "
            f"{poc['total']:>8.2f}s {prod['total']:>9.2f}s {speedup:>6.1f}x | "
            f"{poc['clusters']:>6,}/{prod['clusters']:<7,} "
            f"{poc['matches']:>7,}/{prod['matches']:<8,}"
        )

    print("\nstage breakdown (seconds):")
    print(
        f"  {'dataset':<14} | {'POC idx':>7} {'POC clus':>8} {'POC pair':>8} | "
        f"{'PROD clus*':>10} {'PROD pair':>9}   (*incl. index build)"
    )
    for name, _n, poc, prod in rows:
        print(
            f"  {name:<14} | {poc['index']:>7.2f} {poc['cluster']:>8.2f} "
            f"{poc['pairs']:>8.2f} | {prod['cluster']:>10.2f} {prod['pairs']:>9.2f}"
        )


if __name__ == "__main__":
    main()
