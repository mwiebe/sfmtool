"""Experiment 24 — do long tracks project to 1-D curves in descriptor space?

Claim under test: a track's co-observations (same 3-D point seen across images)
do not form a tight ball in SIFT space; as viewpoint changes continuously the
descriptor moves along a continuous 1-D *curve*. If true, a track's member
descriptors should be:
  - topologically 1-D: their minimum spanning tree is nearly a path (most
    members have MST degree <=2, the diameter path covers ~all members), not a
    star/bush;
  - elongated: the geodesic extent along that path grows ~linearly with track
    length (each new viewpoint extends the arc), not saturating like a ball;
  - curved, not straight: the straight-line (Euclidean) end-to-end distance is
    shorter than the geodesic path length, and global linear (PCA) dimension is
    > 1 even though the local structure is 1-D.

Control: random groups of background (untracked) descriptors of matched size,
which should be high-dimensional bushes (high MST branching, no elongation).

Usage: pixi run -e test python experiments/exp24_track_curves.py [sfmr ...]
"""
from __future__ import annotations
import sys
from collections import defaultdict
import numpy as np
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import minimum_spanning_tree, shortest_path

sys.path.insert(0, "experiments")
from sfm_descriptors import load_descriptor_bank


def track_metrics(X: np.ndarray) -> dict:
    """Curve diagnostics for one group of m descriptors (m, 128) float."""
    m = len(X)
    D = cdist(X, X)  # Euclidean L2, the matcher's space
    mst = minimum_spanning_tree(D)
    W = mst + mst.T
    deg = (W > 0).sum(axis=1).A1 if hasattr((W > 0).sum(axis=1), "A1") else np.asarray((W > 0).sum(axis=1)).ravel()
    # geodesic (weighted) and hop (unweighted) distances along the tree
    geo = shortest_path(W, method="D", directed=False)
    hops = shortest_path((W > 0).astype(float), method="D", directed=False)
    a, b = np.unravel_index(np.argmax(geo), geo.shape)  # diameter endpoints
    geo_diam = geo[a, b]
    path_nodes = int(hops[a, b]) + 1                    # nodes on the diameter path
    edge_w = W.data
    mean_step = float(edge_w.mean())
    eucl_end = float(D[a, b])                           # straight-line end-to-end
    # global linear dimension via PCA participation ratio
    Xc = X - X.mean(0)
    ev = np.linalg.svd(Xc, compute_uv=False) ** 2
    ev = ev[ev > 1e-9]
    pdim = float(ev.sum() ** 2 / (ev ** 2).sum())
    return dict(
        m=m,
        frac_deg_le2=float(np.mean(deg <= 2)),
        branch=int(np.sum(deg >= 3)),
        coverage=path_nodes / m,                        # ~1.0 => MST is a path
        aspect=geo_diam / mean_step,                    # ~m-1 for a clean chain
        straightness=eucl_end / geo_diam if geo_diam else 1.0,  # <1 => curved
        pdim=pdim,
    )


def summarize(label, groups, rng):
    bins = [(3, 4), (5, 7), (8, 12), (13, 20), (21, 999)]
    rows = defaultdict(list)
    ctrl = defaultdict(list)
    bg = groups["__bg__"]
    for pid, X in groups.items():
        if pid == "__bg__":
            continue
        met = track_metrics(X)
        rows[met["m"]].append(met)
        # size-matched background control
        idx = rng.choice(len(bg), size=len(X), replace=False)
        ctrl[met["m"]].append(track_metrics(bg[idx]))

    print(f"\n========== {label} ==========")
    print(f"{'':>9} {'':>5} | {'--- tracks ---':>27} {'straight':>8} {'pdim':>5} "
          f"| {'--- random control ---':>27} {'straight':>8} {'pdim':>5}")
    print(f"{'len bin':>9} {'n':>5} | {'deg<=2':>6} {'branch':>6} {'cover':>5} "
          f"{'aspect':>6} {'straight':>8} {'pdim':>5} | {'deg<=2':>6} {'branch':>6} "
          f"{'cover':>5} {'aspect':>6} {'straight':>8} {'pdim':>5}")
    for lo, hi in bins:
        sel = [met for mm, lst in rows.items() if lo <= mm <= hi for met in lst]
        csel = [met for mm, lst in ctrl.items() if lo <= mm <= hi for met in lst]
        if not sel:
            continue
        def mean(key, s):
            return float(np.mean([x[key] for x in s]))
        name = f"{lo}-{hi if hi < 999 else '+'}"
        def block(s):
            return (f"{mean('frac_deg_le2', s):>6.2f} {mean('branch', s):>6.2f} "
                    f"{mean('coverage', s):>5.2f} {mean('aspect', s):>6.1f} "
                    f"{mean('straightness', s):>8.2f} {mean('pdim', s):>5.1f}")
        print(f"{name:>9} {len(sel):>5} | {block(sel)} | {block(csel)}")


def load_groups(path):
    bank = load_descriptor_bank(path)
    X = bank.descriptors.astype(np.float32)
    pid = bank.point_label
    groups = {}
    for p in np.unique(pid[pid >= 0]):
        rows = np.flatnonzero(pid == p)
        if len(rows) >= 3:
            groups[int(p)] = X[rows]
    groups["__bg__"] = X[pid < 0]
    return groups


if __name__ == "__main__":
    paths = sys.argv[1:] or [
        "/tmp/dino_ab_ws/exhaustive.sfmr",
        "/tmp/seattle_backyard_ab_ws/cluster_recon.sfmr",
    ]
    rng = np.random.default_rng(0)
    for p in paths:
        label = p.split("/")[-2] + ":" + p.split("/")[-1]
        summarize(label, load_groups(p), rng)
