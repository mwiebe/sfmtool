"""Experiment 27 — raw oriented patches as the descriptor; how do they cluster?

Instead of the SIFT gradient histogram, sample a small fronto-parallel patch at
each keypoint: the stored affine-shape matrix is the keypoint's oriented, scaled
local frame, so an 8x8 grid sampled through it (extent f*sigma) is the
orientation-normalized, scale-sized patch. Intensity-normalized (zero-mean,
unit-std) for basic illumination invariance -> a 64-D descriptor.

Compares patch space vs SIFT space on:
  A) intrinsic dimension of tracks vs length (exp24 metric) + random control;
  B) co-obs-vs-background separation AUC on a shared candidate set (sift-space
     k-NN edges), scored by sift L2 vs patch L2.

Usage: pixi run -e test python experiments/exp27_patch_descriptor.py [sfmr ...]
"""
from __future__ import annotations
import sys
import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates
from scipy.stats import rankdata

sys.path.insert(0, "experiments")
from sfm_descriptors import load_descriptor_bank
from sfmtool import KdForest
from sfmtool.sift.file import SiftReader, get_sift_path_for_image

P, F, K = 8, 2.0, 24  # patch grid, extent in sigma, k-NN width


def extract_patches(ws, names):
    g = np.stack(np.meshgrid(np.linspace(-F, F, P), np.linspace(-F, F, P)), -1).reshape(-1, 2)
    out = []
    for name in names:
        im = np.asarray(Image.open(str(ws / name)).convert("L"), dtype=np.float32)
        with SiftReader(get_sift_path_for_image(str(ws / name))) as r:
            pos, aff = r.read_positions_and_shapes()
        pos = np.asarray(pos); aff = np.asarray(aff)        # (K,2), (K,2,2)
        # image coords for every keypoint's PxP grid: center + aff @ g
        disp = np.einsum("kij,gj->kgi", aff, g)             # (K, P*P, 2) (x,y)
        xy = pos[:, None, :] + disp
        rows = xy[..., 1].ravel(); cols = xy[..., 0].ravel()
        vals = map_coordinates(im, [rows, cols], order=1, mode="nearest")
        patches = vals.reshape(len(pos), P * P)
        patches -= patches.mean(1, keepdims=True)
        patches /= patches.std(1, keepdims=True) + 1e-6
        out.append(patches.astype(np.float32))
    return np.concatenate(out, 0)


def pdim(X):
    ev = np.linalg.svd(X - X.mean(0), compute_uv=False) ** 2
    ev = ev[ev > 1e-9]
    return float(ev.sum() ** 2 / (ev ** 2).sum())


def auc(pos, neg):
    r = rankdata(np.concatenate([pos, neg]))
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def run(path):
    label = path.split("/")[-2] + ":" + path.split("/")[-1]
    bank = load_descriptor_bank(path)
    sift = bank.descriptors.astype(np.float32)
    patch = extract_patches(bank.workspace_dir, list(bank.image_names))
    assert len(patch) == len(sift), (len(patch), len(sift))
    img, pid = bank.image_label.astype(np.int64), bank.point_label
    rng = np.random.default_rng(0)
    bg = np.flatnonzero(pid < 0)

    # ---- A) intrinsic dimension vs track length ----
    print(f"\n===== {label}   N={len(sift):,} =====")
    print("  (A) intrinsic dim of tracks vs length  [track | random control]")
    print(f"  {'len bin':>9} {'n':>5} | {'sift pdim':>20} | {'patch pdim':>20}")
    pts, counts = np.unique(pid[pid >= 0], return_counts=True)
    for lo, hi in [(3, 4), (5, 7), (8, 12), (13, 20), (21, 999)]:
        sel = pts[(counts >= lo) & (counts <= hi)]
        if len(sel) == 0:
            continue
        sd, pd, sc, pc = [], [], [], []
        for p in sel:
            rows = np.flatnonzero(pid == p)
            sd.append(pdim(sift[rows])); pd.append(pdim(patch[rows]))
            ridx = rng.choice(bg, size=len(rows), replace=False)
            sc.append(pdim(sift[ridx])); pc.append(pdim(patch[ridx]))
        nm = f"{lo}-{hi if hi < 999 else '+'}"
        print(f"  {nm:>9} {len(sel):>5} | {np.mean(sd):>8.1f} (ctrl {np.mean(sc):>5.1f}) "
              f"| {np.mean(pd):>8.1f} (ctrl {np.mean(pc):>5.1f})")

    # ---- B) separation on sift-space k-NN candidate edges ----
    idx, sdist = KdForest(np.ascontiguousarray(bank.descriptors), preset="accurate").query(
        np.ascontiguousarray(bank.descriptors), k=K + 1)
    idx = idx[:, 1:].astype(np.int64); sdist = sdist[:, 1:]
    seeds = rng.choice(np.flatnonzero(pid >= 0), size=min(4000, int((pid >= 0).sum())), replace=False)
    sp, sn, pp, pn = [], [], [], []
    for i in seeds:
        for c in range(K):
            j = idx[i, c]
            if img[j] == img[i]:
                continue
            co = pid[i] >= 0 and pid[j] == pid[i]
            (sp if co else sn).append(sdist[i, c])
            pdpatch = np.linalg.norm(patch[i] - patch[j])
            (pp if co else pn).append(pdpatch)
    print("  (B) co-obs vs background separation AUC (shared sift-kNN edges):")
    print(f"        sift L2 : {auc([-x for x in sp], [-x for x in sn]):.3f}")
    print(f"        patch L2: {auc([-x for x in pp], [-x for x in pn]):.3f}   "
          f"(co-obs edges {len(sp):,}, background {len(sn):,})")


if __name__ == "__main__":
    for p in sys.argv[1:] or ["/tmp/seattle_backyard_ab_ws/cluster_recon.sfmr"]:
        run(p)
