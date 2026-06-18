"""Visualize long tracks in descriptor space (companion to exp24_track_curves).

For the longest tracks in a reconstruction, PCA-project each track's member
descriptors to 2-D and draw them as a polyline ordered along the track's MST
diameter path. A 1-D curve shows as a smooth, non-self-crossing arc; a ball
shows as a blob with no consistent ordering.

Usage: pixi run -e test python experiments/exp24_viz.py [sfmr] [out.png]
"""
from __future__ import annotations
import sys
import numpy as np
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import minimum_spanning_tree, shortest_path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "experiments")
from sfm_descriptors import load_descriptor_bank


def diameter_order(X):
    D = cdist(X, X)
    W = minimum_spanning_tree(D)
    W = W + W.T
    geo, pred = shortest_path(W, method="D", directed=False, return_predecessors=True)
    a, b = np.unravel_index(np.argmax(geo), geo.shape)
    order, cur = [], b
    while cur != a and cur >= 0:
        order.append(cur)
        cur = pred[a, cur]
    order.append(a)
    return order[::-1]


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dino_ab_ws/exhaustive.sfmr"
    out = sys.argv[2] if len(sys.argv) > 2 else "/tmp/track_curves.png"
    bank = load_descriptor_bank(path)
    X = bank.descriptors.astype(np.float32)
    pid = bank.point_label
    tracks = [(int(p), np.flatnonzero(pid == p)) for p in np.unique(pid[pid >= 0])]
    tracks = [t for t in tracks if len(t[1]) >= 10]
    tracks.sort(key=lambda t: -len(t[1]))
    tracks = tracks[:5]
    bg = np.flatnonzero(pid < 0)
    rng = np.random.default_rng(0)

    def panel(ax, Xt, title):
        order = diameter_order(Xt)
        Xc = Xt - Xt.mean(0)
        _, s, Vt = np.linalg.svd(Xc, full_matrices=False)
        proj = Xc @ Vt[:2].T
        ax.plot(proj[order, 0], proj[order, 1], "-", color="0.6", lw=1, zorder=1)
        ax.scatter(proj[:, 0], proj[:, 1], c=proj[:, 0], cmap="viridis", s=40, zorder=2)
        var2 = (s[:2] ** 2).sum() / (s ** 2).sum()
        ax.set_title(f"{title}: {len(Xt)} pts, PC1+2={var2:.0%} var", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])

    fig, axes = plt.subplots(2, 5, figsize=(18, 7.5))
    for col, (p, rows) in enumerate(tracks):
        panel(axes[0, col], X[rows], f"track {p}")
        ctrl = X[rng.choice(bg, size=len(rows), replace=False)]
        panel(axes[1, col], ctrl, "random bg")
    axes[0, 0].set_ylabel("real tracks", fontsize=12)
    axes[1, 0].set_ylabel("size-matched\nrandom background", fontsize=12)
    fig.suptitle("Descriptor-space structure: longest tracks (top) vs random background groups (bottom)\n"
                 "PCA-2D, points colored by MST-diameter order. Higher PC1+2 variance = lower intrinsic dimension.")
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
