"""Experiment 25 — shared-neighbor (intersection) similarity for clustering.

exp24: raw distances are confounded in 128-D. SNN/Jarvis-Patrick says use the
*overlap* of two descriptors' neighbourhoods instead. Two questions:

  (a) Separation. Among an in-track seed's cross-image k-NN, does the
      neighbourhood-intersection size separate co-observations from background
      better (and more robustly) than raw distance? Reported as AUC for the
      co-obs-vs-background classification, per dataset, distance vs SNN vs
      Jaccard. dino (repetitive, where distance fails most) is the key case.

  (b) Fragment bridging. For same-track member pairs that distance *cannot*
      link (j beyond i's k-NN), is the shared-neighbour overlap still elevated
      above background? If so, SNN edges can stitch the arcs a radius fragments.

Usage: pixi run -e test python experiments/exp25_snn.py [sfmr ...]
"""
from __future__ import annotations
import sys
import numpy as np
from scipy.stats import rankdata

sys.path.insert(0, "experiments")
from sfm_descriptors import load_descriptor_bank
from sfmtool import KdForest

K = 32  # neighbourhood width for SNN


def auc(pos, neg):
    """P(score_pos > score_neg) via rank statistic."""
    pos, neg = np.asarray(pos), np.asarray(neg)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    r = rankdata(np.concatenate([pos, neg]))
    return (r[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def run(path):
    label = path.split("/")[-2] + ":" + path.split("/")[-1]
    bank = load_descriptor_bank(path)
    X = np.ascontiguousarray(bank.descriptors)
    img = bank.image_label.astype(np.int64)
    pid = bank.point_label
    N = len(X)
    idx, dist = KdForest(X, preset="accurate").query(X, k=K + 1)
    idx = idx[:, 1:].astype(np.int64)        # (N, K) neighbour ids, self dropped
    dist = dist[:, 1:]                        # aligned L2 distances
    nbset = [set(row.tolist()) for row in idx]  # neighbour id sets

    def snn(i, j):
        inter = len(nbset[i] & nbset[j])
        return inter, inter / (2 * K - inter) if (2 * K - inter) else 0.0

    # ---- (a) separation on an in-track seed's cross-image neighbours ----
    rng = np.random.default_rng(0)
    in_track = np.flatnonzero(pid >= 0)
    seeds = rng.choice(in_track, size=min(4000, len(in_track)), replace=False)
    d_pos, d_neg, s_pos, s_neg, j_pos, j_neg = [], [], [], [], [], []
    for i in seeds:
        for c in range(K):
            j = idx[i, c]
            if img[j] == img[i]:
                continue
            inter, jac = snn(i, j)
            co = pid[i] >= 0 and pid[j] == pid[i]
            (d_pos if co else d_neg).append(dist[i, c])
            (s_pos if co else s_neg).append(inter)
            (j_pos if co else j_neg).append(jac)

    print(f"\n===== {label}   N={N:,}  co-obs edges={len(s_pos):,}  bg edges={len(s_neg):,} =====")
    print(f"  AUC  (separating co-obs from background, higher=better):")
    print(f"     raw distance : {auc([-x for x in d_pos], [-x for x in d_neg]):.3f}")
    print(f"     SNN |∩|      : {auc(s_pos, s_neg):.3f}")
    print(f"     Jaccard      : {auc(j_pos, j_neg):.3f}")
    print(f"  median shared-neighbour count   co-obs {np.median(s_pos):.0f}  "
          f"vs background {np.median(s_neg):.0f}   (of K={K})")

    # ---- (b) fragment bridging: same-track pairs distance can't link ----
    # near = j in i's k-NN; far = same-track but NOT in i's k-NN.
    near_snn, far_snn = [], []
    pts, counts = np.unique(pid[in_track], return_counts=True)
    long_pts = pts[counts >= 6]
    for p in rng.choice(long_pts, size=min(400, len(long_pts)), replace=False):
        members = np.flatnonzero(pid == p)
        for a_i in range(len(members)):
            for b_i in range(a_i + 1, len(members)):
                i, j = int(members[a_i]), int(members[b_i])
                if img[i] == img[j]:
                    continue
                inter, _ = snn(i, j)
                (near_snn if j in nbset[i] or i in nbset[j] else far_snn).append(inter)
    print(f"  fragment bridging (same-track member pairs):")
    print(f"     near pairs (distance links): median shared {np.median(near_snn):.0f}  "
          f"(n={len(near_snn):,})")
    print(f"     FAR  pairs (distance can't): median shared {np.median(far_snn):.0f}  "
          f"(n={len(far_snn):,})   vs background median {np.median(s_neg):.0f}")
    far = np.asarray(far_snn)
    print(f"     far-pair shared > background-median: {100 * (far > np.median(s_neg)).mean():.0f}%")


if __name__ == "__main__":
    paths = sys.argv[1:] or [
        "/tmp/dino_ab_ws/exhaustive.sfmr",
        "/tmp/seattle_backyard_ab_ws/cluster_recon.sfmr",
    ]
    for p in paths:
        run(p)
