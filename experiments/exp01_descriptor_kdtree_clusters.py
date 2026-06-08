# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 01 — do descriptor-space clusters look like the solve's tracks?

Idea (POC): throw *every* SIFT descriptor from every image into one ANN/KD-tree
and ask whether the structure of descriptor space — nearest neighbours, density
clusters, mutual-NN components — recovers the feature tracks that the SfM solve
actually built.  A "track" here is the set of feature observations the solve tied
to one 3D point.

Four probes, from weakest to strongest claim:

  A. Separation     — are within-track descriptor distances smaller than random
                      cross-image distances?  (Is there any signal at all?)
  B. k-NN recovery  — for a feature in a track, are the *other* members of its
                      track among its descriptor-space nearest neighbours?
  C. DBSCAN vs tracks — density clustering of descriptor space, scored against
                      the track labelling (homogeneity / completeness / ARI).
  D. Mutual-NN comps — cross-image mutual nearest neighbours -> connected
                      components as candidate tracks; how pure / track-like?

Usage:
    pixi run -e experiments python experiments/exp01_descriptor_kdtree_clusters.py \
        seoul_bull_ws/sfmr/*.sfmr
"""

from __future__ import annotations

import argparse
import glob

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.cluster import DBSCAN
from sklearn.metrics import (
    adjusted_rand_score,
    completeness_score,
    homogeneity_score,
)

from knn import knn_indices
from sfm_descriptors import DescriptorBank, load_descriptor_bank

RNG = np.random.default_rng(42)

# SIFT descriptors are 128-D.  A literal KD-tree is pathological at this
# dimensionality (it degrades to worse-than-linear), so the "kd-tree" in the
# experiment's name is really "exact nearest neighbours" — computed via a
# chunked BLAS GEMM (see knn.py), which scales to the ~340k-descriptor datasets.
# A production system would swap in a real ANN index (faiss / HNSW).

# DBSCAN over *all* descriptors doesn't scale; above this many we sample (keeping
# every in-track descriptor) so the density clustering probe stays tractable.
# Large eps on a big subsample blows up the neighbour graph (OOM), so we also
# cap the eps sweep below the always-collapses-to-mega-clusters regime.
DBSCAN_MAX_N = 80_000
DBSCAN_EPS_PERCENTILES = (25, 50, 75)


def _header(title: str) -> None:
    print(f"\n{'=' * 70}\n{title}\n{'=' * 70}")


def probe_a_separation(bank: DescriptorBank, n_pairs: int = 30_000) -> np.ndarray:
    """Compare within-track vs random cross-image descriptor distances.

    Returns the within-track distance samples so probe C can sweep eps over
    their percentiles.
    """
    _header("A. Separation: within-track vs random cross-image distances")
    X = bank.descriptors.astype(np.float32)
    pid = bank.point_label
    img = bank.image_label

    # Within-track pairs: sample pairs of observations sharing a point id.
    in_track = np.flatnonzero(pid >= 0)
    order = np.argsort(pid[in_track], kind="stable")
    sorted_rows = in_track[order]
    sorted_pid = pid[sorted_rows]
    # group boundaries
    bounds = np.flatnonzero(np.diff(sorted_pid)) + 1
    groups = np.split(sorted_rows, bounds)

    within = []
    for g in groups:
        if len(g) < 2:
            continue
        # all-pairs within a small track is cheap
        for a in range(len(g)):
            for b in range(a + 1, len(g)):
                within.append((g[a], g[b]))
    within = np.asarray(within)
    if len(within) > n_pairs:
        within = within[RNG.choice(len(within), n_pairs, replace=False)]
    wd = np.linalg.norm(X[within[:, 0]] - X[within[:, 1]], axis=1)

    # Random cross-image pairs.
    ri = RNG.integers(0, bank.n, size=n_pairs * 2)
    rj = RNG.integers(0, bank.n, size=n_pairs * 2)
    keep = img[ri] != img[rj]
    ri, rj = ri[keep][:n_pairs], rj[keep][:n_pairs]
    rd = np.linalg.norm(X[ri] - X[rj], axis=1)

    def stats(d):
        return (
            f"n={len(d):>6}  min={d.min():6.1f}  "
            f"p05={np.percentile(d, 5):6.1f}  median={np.median(d):6.1f}  "
            f"p95={np.percentile(d, 95):6.1f}  max={d.max():6.1f}"
        )

    print(f"  within-track : {stats(wd)}")
    print(f"  random pairs : {stats(rd)}")
    # Crude separability: how often is a within-track distance below the random p05?
    thr = np.percentile(rd, 5)
    frac_below = float(np.mean(wd < thr))
    print(f"  within-track distances below random p05 ({thr:.1f}): {frac_below:.1%}")
    return wd


def probe_b_knn_recovery(bank: DescriptorBank, k: int = 5) -> None:
    """For in-track features, are co-track members among their k-NN?"""
    _header(f"B. k-NN track recovery (k={k}, cross-image neighbours only)")
    X = bank.descriptors.astype(np.float32)
    pid = bank.point_label
    img = bank.image_label

    # Over-fetch so we can drop self + same-image neighbours and still keep k.
    n_fetch = k + 40
    q = np.flatnonzero(pid >= 0)
    idx = knn_indices(X, X[q], min(n_fetch, bank.n))

    precisions = []
    recalls = []
    track_sizes = {p: c for p, c in zip(*np.unique(pid[q], return_counts=True))}
    for row, qi in enumerate(q):
        neigh = idx[row]
        neigh = neigh[neigh != qi]
        neigh = neigh[img[neigh] != img[qi]]  # matches live in *other* images
        neigh = neigh[:k]
        if len(neigh) == 0:
            continue
        same = pid[neigh] == pid[qi]
        precisions.append(float(np.mean(same)))
        # recall vs the other members of this track
        others = track_sizes[pid[qi]] - 1
        if others > 0:
            recalls.append(float(np.sum(same)) / others)

    print(f"  queried in-track features : {len(q)}")
    print(f"  mean precision@{k}        : {np.mean(precisions):.3f}")
    print(f"  mean recall@{k}           : {np.mean(recalls):.3f}")
    # fraction of features whose #1 cross-image neighbour is a true co-track member
    hit1 = 0
    for row, qi in enumerate(q):
        neigh = idx[row]
        neigh = neigh[neigh != qi]
        neigh = neigh[img[neigh] != img[qi]]
        if len(neigh) and pid[neigh[0]] == pid[qi]:
            hit1 += 1
    print(f"  top-1 cross-image neighbour is co-track: {hit1 / len(q):.1%}")


def probe_c_dbscan(
    bank: DescriptorBank, within: np.ndarray, eps_override: float | None = None
) -> None:
    """DBSCAN over descriptor space, swept across eps, scored against tracks.

    eps is the make-or-break knob: too large and every dense blob merges into a
    few mega-clusters; too small and real tracks fragment.  We sweep the lower
    percentiles of the within-track distance distribution.
    """
    _header("C. DBSCAN vs tracks (eps swept over within-track percentiles)")
    n_solve_tracks = len(np.unique(bank.point_label[bank.point_label >= 0]))
    print(f"  solve tracks: {n_solve_tracks}")

    # DBSCAN doesn't scale to hundreds of thousands of points; subsample the
    # corpus but keep every in-track descriptor so the scoring is unaffected.
    if bank.n > DBSCAN_MAX_N:
        in_track = np.flatnonzero(bank.point_label >= 0)
        rest = np.flatnonzero(bank.point_label < 0)
        n_fill = max(0, DBSCAN_MAX_N - len(in_track))
        fill = RNG.choice(rest, size=min(n_fill, len(rest)), replace=False)
        sel = np.sort(np.concatenate([in_track, fill]))
        print(
            f"  (subsampled {len(sel)} of {bank.n} descriptors; all "
            f"{len(in_track)} in-track kept)"
        )
    else:
        sel = np.arange(bank.n)

    X = bank.descriptors[sel].astype(np.float32)
    point_label = bank.point_label[sel]

    if eps_override is not None:
        eps_list = [(None, eps_override)]
    else:
        eps_list = [
            (p, float(np.percentile(within, p))) for p in DBSCAN_EPS_PERCENTILES
        ]

    mask = point_label >= 0
    gt = point_label[mask]
    print(f"  scored on the {mask.sum()} in-track descriptors:\n")
    print(
        f"  {'eps':>22}  {'clusters':>8}  {'noise':>6}  "
        f"{'homog':>6}  {'compl':>6}  {'ARI':>6}"
    )
    for pct, eps in eps_list:
        labels = DBSCAN(eps=eps, min_samples=2, n_jobs=-1).fit_predict(X)
        n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise = int(np.sum(labels == -1))
        pred = labels[mask]
        tag = f"p{pct}={eps:6.1f}" if pct is not None else f"{eps:6.1f}"
        print(
            f"  {tag:>22}  {n_clusters:>8}  {n_noise:>6}  "
            f"{homogeneity_score(gt, pred):>6.3f}  "
            f"{completeness_score(gt, pred):>6.3f}  "
            f"{adjusted_rand_score(gt, pred):>6.3f}"
        )


def probe_d_mutual_nn(bank: DescriptorBank) -> None:
    """Cross-image mutual-NN graph -> connected components as candidate tracks."""
    _header("D. Mutual-NN connected components as candidate tracks")
    X = bank.descriptors.astype(np.float32)
    img = bank.image_label
    pid = bank.point_label

    # Nearest cross-image neighbour for every descriptor.
    n_fetch = 30
    idx = knn_indices(X, X, min(n_fetch, bank.n))

    # First cross-image neighbour per descriptor (vectorised): mask self and
    # same-image columns, take the first surviving column of each row.
    rng = np.arange(bank.n)
    self_or_same = (idx == rng[:, None]) | (img[idx] == img[:, None])
    valid = ~self_or_same
    first_col = np.argmax(valid, axis=1)
    has_cross = valid[rng, first_col]
    nearest_cross = np.where(has_cross, idx[rng, first_col], -1)

    # Mutual edges only: i->j and j->i, kept once with i < j.
    j = nearest_cross
    has_j = j >= 0
    back = np.full(bank.n, -1, dtype=np.int64)
    back[has_j] = nearest_cross[j[has_j]]
    mutual = has_j & (back == rng) & (rng < j)
    rows = rng[mutual]
    cols = j[mutual]
    n_edges = len(rows)
    data = np.ones(n_edges, dtype=np.int8)
    graph = csr_matrix((data, (rows, cols)), shape=(bank.n, bank.n))
    n_comp, comp = connected_components(graph, directed=False)

    # Keep only non-trivial components (size >= 2).
    sizes = np.bincount(comp)
    multi = np.flatnonzero(sizes >= 2)
    cand_tracks = len(multi)
    print(f"  mutual-NN edges            : {n_edges}")
    print(f"  candidate tracks (size>=2) : {cand_tracks}")
    print(f"  solve tracks               : {len(np.unique(pid[pid >= 0]))}")

    # Group descriptors by component once (avoids O(n) scan per component).
    multi_set = set(multi.tolist())
    order = np.argsort(comp, kind="stable")
    comp_sorted = comp[order]
    bounds = np.flatnonzero(np.diff(comp_sorted)) + 1
    groups = np.split(order, bounds)

    purities = []
    img_unique = []
    covered_in_track = 0
    total_in_track_members = 0
    for members in groups:
        if comp[members[0]] not in multi_set:
            continue
        labs = pid[members]
        in_t = labs[labs >= 0]
        total_in_track_members += len(in_t)
        if len(in_t):
            _, counts = np.unique(in_t, return_counts=True)
            purities.append(counts.max() / len(members))
            covered_in_track += counts.max()
        # real tracks never have two features from one image
        img_unique.append(len(np.unique(img[members])) == len(members))

    if purities:
        print(f"  mean component purity      : {np.mean(purities):.3f}")
    print(f"  components with unique-image members: {np.mean(img_unique):.1%}")
    if total_in_track_members:
        print(
            f"  dominant-track members / in-track members: "
            f"{covered_in_track / total_in_track_members:.1%}"
        )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("sfmr", help="path (or glob) to a .sfmr reconstruction")
    ap.add_argument("-k", type=int, default=5, help="k for k-NN recovery probe")
    ap.add_argument(
        "--eps",
        type=float,
        default=None,
        help="DBSCAN eps (default: derived from probe A)",
    )
    args = ap.parse_args()

    matches = sorted(glob.glob(args.sfmr))
    if not matches:
        raise SystemExit(f"no .sfmr matched: {args.sfmr}")
    path = matches[0]

    print(f"Loading descriptors + tracks from: {path}")
    bank = load_descriptor_bank(path)
    n_in = int(np.sum(bank.point_label >= 0))
    print(f"  images           : {len(bank.image_names)}")
    print(f"  total descriptors: {bank.n}")
    print(f"  in solve tracks  : {n_in} ({n_in / bank.n:.1%})")
    print(
        f"  solve tracks     : {len(np.unique(bank.point_label[bank.point_label >= 0]))}"
    )

    within = probe_a_separation(bank)
    probe_b_knn_recovery(bank, k=args.k)
    # Probe D (the most informative) before C: DBSCAN is the memory-heavy one.
    probe_d_mutual_nn(bank)
    probe_c_dbscan(bank, within, eps_override=args.eps)
    print()


if __name__ == "__main__":
    main()
