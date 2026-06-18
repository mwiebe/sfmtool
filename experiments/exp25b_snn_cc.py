"""Experiment 25b — does SNN-graph clustering find LARGER clusters than the floor?

exp25 showed shared-neighbour overlap separates co-obs from background robustly
(better than distance on dino) but does NOT pairwise-bridge fragments (far
same-track pairs share ~background). Open question: can *transitive* SNN
connectivity still stitch the consecutive arcs into fewer, larger clusters —
and how much does it over-merge distinct points doing so?

Build the shared-neighbour graph over in-track descriptors (edge iff cross-image
k-NN with shared >= tau), take connected components, and compare to the
background-floor clusters on the same reconstruction:
  - fragmentation: distinct clusters a true track's members split across (lower
    = larger/less fragmented; the floor sits ~1.6-2.7).
  - over-merge: clusters mixing >=2 distinct true points (purity).

Usage: pixi run -e test python experiments/exp25b_snn_cc.py [sfmr ...]
"""
from __future__ import annotations
import sys
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, "experiments")
from sfm_descriptors import load_descriptor_bank
from sfmtool import KdForest, background_floor_clusters

K = 32


def frag_purity(comp_of, pid_sub):
    """comp_of: component id per node; pid_sub: true point id per node (>=0)."""
    # fragmentation: distinct components per true point
    frags, mixes = [], []
    for p in np.unique(pid_sub):
        c = np.unique(comp_of[pid_sub == p])
        frags.append(len(c))
    # over-merge: distinct true points per component
    for c in np.unique(comp_of):
        pts = np.unique(pid_sub[comp_of == c])
        mixes.append(len(pts))
    return np.mean(frags), np.mean(np.asarray(mixes) > 1), np.max(mixes)


def run(path):
    label = path.split("/")[-2] + ":" + path.split("/")[-1]
    bank = load_descriptor_bank(path)
    X = np.ascontiguousarray(bank.descriptors)
    img = bank.image_label.astype(np.int64)
    pid = bank.point_label
    starts = bank.image_starts.astype(np.uint32)
    idx, _ = KdForest(X, preset="accurate").query(X, k=K + 1)
    idx = idx[:, 1:].astype(np.int64)
    nbset = [set(r.tolist()) for r in idx]
    print(f"\n===== {label}  N={len(X):,} =====")

    # ---- floor baseline ----
    cs, mi, mf = background_floor_clusters(X, starts, d=10, alpha=0.8)
    rows = starts[mi.astype(np.int64)] + mf.astype(np.int64)
    cl_of_row = np.full(len(X), -1, np.int64)
    for c in range(len(cs) - 1):
        cl_of_row[rows[cs[c]:cs[c + 1]]] = c
    it = (pid >= 0) & (cl_of_row >= 0)
    f, ov, mx = frag_purity(cl_of_row[it], pid[it])
    print(f"  floor (d=10,a=0.8): clusters={len(cs)-1:,}  frag/track={f:.2f}  "
          f"over-merged clusters={ov:.1%}  worst-mix={mx}")

    # ---- SNN connected components over in-track nodes, sweeping tau ----
    in_track = np.flatnonzero(pid >= 0)
    pos = {int(g): k for k, g in enumerate(in_track)}  # global row -> compact id
    pid_sub = pid[in_track]
    for tau in (3, 6, 10, 14):
        ea, eb = [], []
        for k, i in enumerate(in_track):
            for j in idx[i]:
                if j in pos and img[j] != img[i] and len(nbset[i] & nbset[int(j)]) >= tau:
                    ea.append(k); eb.append(pos[int(j)])
        if not ea:
            print(f"  SNN tau={tau:>2}: no edges"); continue
        g = coo_matrix((np.ones(len(ea)), (ea, eb)), shape=(len(in_track),) * 2)
        n_c, comp = connected_components(g + g.T, directed=False)
        f, ov, mx = frag_purity(comp, pid_sub)
        print(f"  SNN tau={tau:>2}: components={n_c:,}  frag/track={f:.2f}  "
              f"over-merged clusters={ov:.1%}  worst-mix={mx}  edges={len(ea):,}")


if __name__ == "__main__":
    for p in sys.argv[1:] or ["/tmp/dino_ab_ws/exhaustive.sfmr",
                              "/tmp/seattle_backyard_ab_ws/cluster_recon.sfmr"]:
        run(p)
