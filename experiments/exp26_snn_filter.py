"""Experiment 26 — SNN as a pre-verification quality filter on floor matches.

The floor matcher emits candidate matches that then go through geometric
verification (RANSAC) and into the solve. Most background candidates are
rejected by verification; if a cheap shared-neighbour test can drop them first,
it offloads RANSAC and may clean the solve's input.

For each floor candidate edge we compute SNN = |N_i ∩ N_j| (shared neighbours,
from an all-descriptor k-NN) and label it co-observation / background via the
reference solve's tracks. Sweeping a threshold tau, we report how much
background it removes (speed/cleanliness) against how many co-observations it
keeps (quality), plus the image-pair reduction (verification runs per pair).

Caveat: "background" = endpoints not in the same reference track. The reference
solve is lean, so some background is actually correct-but-untracked (exp19) —
making background-removed an upper bound and co-obs-kept a lower bound.

Usage: pixi run -e test python experiments/exp26_snn_filter.py [sfmr ...]
"""
from __future__ import annotations
import sys
import numpy as np

sys.path.insert(0, "experiments")
from sfm_descriptors import load_descriptor_bank
from sfmtool import KdForest
from sfmtool._sfmtool import background_floor_clusters, clusters_to_pair_matches

K = 32


def run(path):
    label = path.split("/")[-2] + ":" + path.split("/")[-1]
    bank = load_descriptor_bank(path)
    X = np.ascontiguousarray(bank.descriptors)
    starts = bank.image_starts.astype(np.uint32)
    pid = bank.point_label

    idx, _ = KdForest(X, preset="accurate").query(X, k=K + 1)
    idx = idx[:, 1:].astype(np.int64)
    nb = [frozenset(r.tolist()) for r in idx]

    # floor candidate matches (pre-verification)
    cs, mi, mf = background_floor_clusters(X, starts, d=10, alpha=0.8, min_size=2)
    pairs, counts, feat_pairs, _ = clusters_to_pair_matches(cs, mi, mf, X, starts)
    pairs, counts, feat_pairs = np.asarray(pairs), np.asarray(counts), np.asarray(feat_pairs)
    imgA = np.repeat(pairs[:, 0], counts).astype(np.int64)
    imgB = np.repeat(pairs[:, 1], counts).astype(np.int64)
    gA = starts[imgA] + feat_pairs[:, 0]
    gB = starts[imgB] + feat_pairs[:, 1]
    pair_id = np.repeat(np.arange(len(pairs)), counts)

    coobs = (pid[gA] >= 0) & (pid[gA] == pid[gB])
    snn = np.fromiter((len(nb[a] & nb[b]) for a, b in zip(gA.tolist(), gB.tolist())),
                      dtype=np.int32, count=len(gA))

    M = len(snn); nco = int(coobs.sum()); nbg = M - nco
    print(f"\n===== {label} =====")
    print(f"  floor candidates: {M:,} edges over {len(pairs):,} image pairs "
          f"| co-obs {nco:,} ({100*nco/M:.0f}%)  background {nbg:,} ({100*nbg/M:.0f}%)")
    print(f"  {'tau':>4} {'edges kept':>11} {'co-obs kept':>12} {'bg removed':>11} "
          f"{'pairs kept':>11}")
    for tau in (1, 2, 4, 6, 8, 10):
        keep = snn >= tau
        cok = 100 * (keep & coobs).sum() / nco
        bgr = 100 * (~keep & ~coobs).sum() / nbg
        pk = len(np.unique(pair_id[keep]))
        print(f"  {tau:>4} {keep.sum():>11,} {cok:>11.1f}% {bgr:>10.1f}% "
              f"{pk:>7,} ({100*pk/len(pairs):>3.0f}%)")


if __name__ == "__main__":
    for p in sys.argv[1:] or ["/tmp/dino_ab_ws/exhaustive.sfmr",
                              "/tmp/seattle_backyard_ab_ws/cluster_recon.sfmr"]:
        run(p)
