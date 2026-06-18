"""Plot the backing data for exp27 (patch vs SIFT clustering).

Produces one figure:
  row 0: intrinsic dimension vs track length (track & random control,
         SIFT & patch) for dino and seattle -- the core claim.
  row 1: ROC for co-obs-vs-background separation, SIFT vs patch.
  row 2: one real long dino track drawn in SIFT space vs patch space
         (PCA-2D, colored by the SAME SIFT-chain order).
  row 3: that track's actual 8x8 patches in chain order -- appearance morphing.

Usage: pixi run -e test python experiments/exp27_plot.py
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
from sfmtool import KdForest
from exp27_patch_descriptor import extract_patches, pdim, P, K

BINS = [(3, 4), (5, 7), (8, 12), (13, 20), (21, 999)]
CENTERS = [3.5, 6, 10, 16, 28]


def diameter_info(Xt):
    """Return (diameter-path node order, per-node geodesic coord from one end)."""
    W = minimum_spanning_tree(cdist(Xt, Xt))
    W = W + W.T
    geo, pred = shortest_path(W, directed=False, return_predecessors=True)
    a, b = np.unravel_index(np.argmax(geo), geo.shape)
    order, cur = [], b
    while cur != a and cur >= 0:
        order.append(cur); cur = pred[a, cur]
    order.append(a)
    return order[::-1], geo[a]  # coord covers ALL nodes


def gather(path):
    bank = load_descriptor_bank(path)
    sift = bank.descriptors.astype(np.float32)
    patch = extract_patches(bank.workspace_dir, list(bank.image_names))
    img, pid = bank.image_label.astype(np.int64), bank.point_label
    rng = np.random.default_rng(0)
    bg = np.flatnonzero(pid < 0)
    pts, counts = np.unique(pid[pid >= 0], return_counts=True)

    dim = {k: [] for k in ("st", "sc", "pt", "pc")}  # sift/patch track/ctrl
    for lo, hi in BINS:
        sel = pts[(counts >= lo) & (counts <= hi)]
        st = sc = pt = pc = np.nan
        if len(sel):
            a = b = c = d = []
            a = [pdim(sift[np.flatnonzero(pid == p)]) for p in sel]
            b = [pdim(patch[np.flatnonzero(pid == p)]) for p in sel]
            ctl = [rng.choice(bg, size=int(n), replace=False)
                   for n in counts[(counts >= lo) & (counts <= hi)]]
            c = [pdim(sift[r]) for r in ctl]
            d = [pdim(patch[r]) for r in ctl]
            st, pt, sc, pc = np.mean(a), np.mean(b), np.mean(c), np.mean(d)
        dim["st"].append(st); dim["pt"].append(pt); dim["sc"].append(sc); dim["pc"].append(pc)

    # ROC edges
    idx, sd = KdForest(np.ascontiguousarray(bank.descriptors), preset="accurate").query(
        np.ascontiguousarray(bank.descriptors), k=K + 1)
    idx = idx[:, 1:].astype(np.int64); sd = sd[:, 1:]
    seeds = rng.choice(np.flatnonzero(pid >= 0), size=min(2500, int((pid >= 0).sum())), replace=False)
    co, sdist, pdist = [], [], []
    for i in seeds:
        for c in range(K):
            j = idx[i, c]
            if img[j] == img[i]:
                continue
            co.append(pid[i] == pid[j])
            sdist.append(sd[i, c])
            pdist.append(np.linalg.norm(patch[i] - patch[j]))
    return bank, sift, patch, pid, dim, (np.array(co), np.array(sdist), np.array(pdist))


def roc(co, dist):
    s = -dist  # smaller distance => more co-obs
    o = np.argsort(-s)
    co = co[o]
    tp = np.cumsum(co) / co.sum()
    fp = np.cumsum(~co) / (~co).sum()
    auc = np.sum((fp[1:] - fp[:-1]) * (tp[1:] + tp[:-1]) / 2)
    return fp, tp, auc


def main():
    data = {n: gather(p) for n, p in
            [("dino", "/tmp/dino_ab_ws/exhaustive.sfmr"),
             ("seattle", "/tmp/seattle_backyard_ab_ws/cluster_recon.sfmr")]}

    fig = plt.figure(figsize=(15, 17))
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1.1, 0.7], hspace=0.32, wspace=0.22)

    # row 0: intrinsic dim vs length
    for col, name in enumerate(("dino", "seattle")):
        ax = fig.add_subplot(gs[0, col])
        dim = data[name][4]
        ax.plot(CENTERS, dim["st"], "o-", color="C0", label="SIFT — track")
        ax.plot(CENTERS, dim["sc"], "o--", color="C0", alpha=0.45, label="SIFT — random ctrl")
        ax.plot(CENTERS, dim["pt"], "s-", color="C3", label="patch — track")
        ax.plot(CENTERS, dim["pc"], "s--", color="C3", alpha=0.45, label="patch — random ctrl")
        ax.set_title(f"{name}: intrinsic dimension vs track length")
        ax.set_xlabel("track length (views)"); ax.set_ylabel("PCA participation dim")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # row 1: ROC
    for col, name in enumerate(("dino", "seattle")):
        ax = fig.add_subplot(gs[1, col])
        co, sdist, pdist = data[name][5]
        for dist, c, lab in ((sdist, "C0", "SIFT L2"), (pdist, "C3", "patch L2")):
            fp, tp, a = roc(co, dist)
            ax.plot(fp, tp, color=c, label=f"{lab}  AUC={a:.3f}")
        ax.plot([0, 1], [0, 1], "k:", alpha=0.4)
        ax.set_title(f"{name}: co-obs vs background separation")
        ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate")
        ax.legend(fontsize=9); ax.grid(alpha=0.3)

    # rows 2-3: one long dino track in both spaces + its patches
    bank, sift, patch, pid, _, _ = data["dino"]
    pts, counts = np.unique(pid[pid >= 0], return_counts=True)
    cand = pts[(counts >= 18) & (counts <= 30)]
    p = cand[np.argmax([pdim(patch[np.flatnonzero(pid == q)]) for q in cand])]
    rows = np.flatnonzero(pid == p)
    order, coord = diameter_info(sift[rows])  # SIFT-chain order + coord, used for BOTH

    for col, (Xspace, c, tag) in enumerate(
            ((sift[rows], "C0", "SIFT"), (patch[rows], "C3", "patch"))):
        ax = fig.add_subplot(gs[2, col])
        Xc = Xspace - Xspace.mean(0)
        _, s, Vt = np.linalg.svd(Xc, full_matrices=False)
        proj = Xc @ Vt[:2].T
        ax.plot(proj[order, 0], proj[order, 1], "-", color="0.7", lw=1, zorder=1)
        ax.scatter(proj[:, 0], proj[:, 1], c=coord, cmap="viridis", s=55, zorder=2)
        var2 = (s[:2] ** 2).sum() / (s ** 2).sum()
        ax.set_title(f"track {p} ({len(rows)} views) in {tag} space — "
                     f"pdim={pdim(Xspace):.1f}, PC1+2={var2:.0%} var")
        ax.set_xticks([]); ax.set_yticks([])

    # patches as individual SQUARE cells, in chain order
    cap = fig.add_subplot(gs[3, :]); cap.axis("off")
    cap.set_title(f"track {p}: its actual 8x8 oriented patches, in chain order "
                  "(appearance morphs across viewpoints)", fontsize=11)
    po = patch[rows][np.argsort(coord)]
    n = len(rows); ncol = min(n, 14); nrow = int(np.ceil(n / ncol))
    sub = gs[3, :].subgridspec(nrow, ncol, wspace=0.08, hspace=0.08)
    for i in range(n):
        a = fig.add_subplot(sub[i // ncol, i % ncol])
        a.imshow(po[i].reshape(P, P), cmap="gray")  # imshow default aspect='equal' => square
        a.set_xticks([]); a.set_yticks([])

    fig.suptitle("Raw oriented patches vs SIFT as a clustering descriptor\n"
                 "SIFT's invariance compresses a track's cross-view variation into a low-D curve; "
                 "raw patches let it spread into a high-D blob.", fontsize=13)
    out = "/tmp/patch_vs_sift.png"
    fig.savefig(out, dpi=110, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
