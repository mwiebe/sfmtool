"""Show a track's oriented patches in COLOR (companion to exp27).

Samples RGB patches through each keypoint's affine-shape frame (same frame as
the descriptor, rendered at finer resolution for display) and lays out several
long dino tracks, one row each, in SIFT-chain order. Shows the cross-view
appearance morphing in colour.

Usage: pixi run -e test python experiments/exp27_color_patches.py
"""
from __future__ import annotations
import sys
import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import minimum_spanning_tree, shortest_path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "experiments")
from sfm_descriptors import load_descriptor_bank
from sfmtool.sift.file import SiftReader, get_sift_path_for_image

F, PV = 2.0, 32  # extent in sigma; PV = display sampling resolution


def chain_order(Xt):
    W = minimum_spanning_tree(cdist(Xt, Xt)); W = W + W.T
    geo, pred = shortest_path(W, directed=False, return_predecessors=True)
    a, b = np.unravel_index(np.argmax(geo), geo.shape)
    order, cur = [], b
    while cur != a and cur >= 0:
        order.append(cur); cur = pred[a, cur]
    order.append(a)
    return order[::-1]


def main():
    path = "/tmp/dino_ab_ws/exhaustive.sfmr"
    bank = load_descriptor_bank(path)
    sift = bank.descriptors.astype(np.float32)
    img_label = bank.image_label.astype(np.int64)
    feat_label = bank.feature_label.astype(np.int64)
    pid = bank.point_label
    ws, names = bank.workspace_dir, list(bank.image_names)

    pts, counts = np.unique(pid[pid >= 0], return_counts=True)
    longest = pts[np.argsort(-counts)]
    pick = [6019] if 6019 in set(pts.tolist()) else []
    for q in longest:
        if len(pick) >= 5:
            break
        if int(q) not in pick and counts[pts == q][0] <= 34:
            pick.append(int(q))

    g = np.stack(np.meshgrid(np.linspace(-F, F, PV), np.linspace(-F, F, PV)), -1).reshape(-1, 2)
    img_cache, kp_cache = {}, {}

    def patch_rgb(global_row):
        im_i = img_label[global_row]
        if im_i not in img_cache:
            img_cache[im_i] = np.asarray(Image.open(str(ws / names[im_i])).convert("RGB"), np.float32)
            with SiftReader(get_sift_path_for_image(str(ws / names[im_i]))) as r:
                kp_cache[im_i] = r.read_positions_and_shapes()
        im = img_cache[im_i]
        pos, aff = kp_cache[im_i]
        f = feat_label[global_row]
        xy = np.asarray(pos[f]) + (np.asarray(aff[f]) @ g.T).T  # (PV*PV, 2)
        chans = [map_coordinates(im[..., c], [xy[:, 1], xy[:, 0]], order=1, mode="nearest")
                 for c in range(3)]
        patch = np.stack(chans, -1).reshape(PV, PV, 3)
        patch = (patch - patch.min()) / (np.ptp(patch) + 1e-6)  # per-patch contrast stretch
        return patch

    ncol = 16
    fig, axes = plt.subplots(len(pick), ncol, figsize=(ncol, len(pick) + 0.6),
                             squeeze=False)
    for r, p in enumerate(pick):
        rows = np.flatnonzero(pid == p)
        order = chain_order(sift[rows])
        rows = rows[order][:ncol]
        for c in range(ncol):
            ax = axes[r][c]; ax.set_xticks([]); ax.set_yticks([])
            if c < len(rows):
                ax.imshow(patch_rgb(rows[c]))
            else:
                ax.axis("off")
        axes[r][0].set_ylabel(f"track {p}\n{len(np.flatnonzero(pid == p))} views",
                              fontsize=8, rotation=0, ha="right", va="center")
    fig.suptitle("Oriented colour patches of long dino tracks, in chain order "
                 "(one row per track) — the same surface point across viewpoints", fontsize=11)
    fig.tight_layout()
    out = "/tmp/track_patches_color.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
