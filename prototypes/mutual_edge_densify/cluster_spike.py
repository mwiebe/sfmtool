"""Spike (cluster granularity): mutual-kNN *clusters* as track candidates.

Companion to spike.py, which tested mutual edges individually (two-view). This
evaluates the natural unit instead: a mutual-kNN cluster -- the connected
component of mutual edges over (image, feature) nodes -- is a putative
multi-view track. Each cluster is triangulated as a whole (N-view), outlier
members are trimmed, and the surviving cluster is reconciled against the
existing reconstruction's tracks.

Reports actual candidate *tracks* (new / extended / merged), the over-merge
diagnostics edge-level can't see (cluster-size distribution, clusters with two
features in the same image), and cost.

Run: pixi run -e test python prototypes/mutual_edge_densify/cluster_spike.py
"""
from __future__ import annotations
import sys, time
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from sfmtool._sfmtool import (
    SfmrReconstruction, mutual_knn_matches, triangulate_batch, read_sift)

REPROJ_PX = 4.0      # max per-member reprojection error to keep a member
MIN_ANGLE_DEG = 1.5  # min triangulation spread to accept a cluster (proxy)
K = 20


def quat_wxyz_to_R(q):  # world->cam, COLMAP convention
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ])


def triangulate_and_reproject(m_img, m_feat, offsets, rays_world, centers, kp, R, t, cams, N):
    """N-view triangulate each group; return point, per-member reproj err + depth."""
    dirs = np.stack([rays_world[i][f] for i, f in zip(m_img, m_feat)])
    cen = centers[m_img]
    tri = triangulate_batch(dirs, cen, offsets)
    pts = np.asarray(tri["points"]); cond = np.asarray(tri["condition_number"])
    sizes = np.diff(offsets)
    Xmem = np.repeat(pts, sizes, axis=0)
    px = np.stack([kp[i][f] for i, f in zip(m_img, m_feat)])
    err = np.full(len(m_img), np.inf); depth = np.full(len(m_img), -1.0)
    for n in range(N):
        mask = m_img == n
        if mask.any():
            xc = (R[n] @ Xmem[mask].T).T + t[n]
            depth[mask] = xc[:, 2]
            pr = np.asarray(cams[n].ray_to_pixel_batch(xc))
            err[mask] = np.linalg.norm(pr - px[mask], axis=1)
    # triangulation spread proxy: 2 * half-angle from the cluster mean ray
    mean = np.add.reduceat(dirs, offsets[:-1].astype(np.intp), axis=0)
    mean /= np.linalg.norm(mean, axis=1, keepdims=True) + 1e-12
    cosm = np.sum(dirs * np.repeat(mean, sizes, axis=0), axis=1)
    spread_deg = 2 * np.degrees(np.arccos(np.clip(
        np.minimum.reduceat(cosm, offsets[:-1].astype(np.intp)), -1, 1)))
    return pts, cond, err, depth, spread_deg


def run(name, recon_path, features_dir, triangle_mins=(0,)):
    print(f"\n========== {name} ==========", flush=True)
    recon = SfmrReconstruction.load(str(recon_path))
    names = recon.image_names
    N = len(names)
    R = np.stack([quat_wxyz_to_R(q) for q in recon.quaternions_wxyz])
    t = np.asarray(recon.translations)
    centers = np.einsum("nij,nj->ni", R.transpose(0, 2, 1), -t)
    cam_idx = recon.camera_indexes
    cams = [recon.cameras[int(c)] for c in cam_idx]
    positions = recon.positions
    ti, tf, tp = recon.track_image_indexes, recon.track_feature_indexes, recon.track_point_ids

    kp, descs, starts = [], [], [0]
    for nm in names:
        d = read_sift(str(features_dir / f"{Path(nm).name}.sift"))
        kp.append(np.asarray(d["positions_xy"], dtype=np.float64))
        descs.append(np.asarray(d["descriptors"], dtype=np.uint8))
        starts.append(starts[-1] + len(descs[-1]))
    corpus = np.concatenate(descs, axis=0)
    image_starts = np.asarray(starts, dtype=np.uint32)
    rays_world = [(R[n].T @ np.asarray(cams[n].pixel_to_ray_batch(kp[n])).T).T for n in range(N)]

    # pose self-check
    errs = []
    for n in range(N):
        m = ti == n
        if m.any():
            xc = (R[n] @ positions[tp[m]].T).T + t[n]
            e = np.linalg.norm(np.asarray(cams[n].ray_to_pixel_batch(xc)) - kp[n][tf[m]], axis=1)
            errs.append(e[np.isfinite(e)])
    med = float(np.median(np.concatenate(errs)))
    print(f"pose self-check: median reproj {med:.3f}px -> {'OK' if med < 2 else 'BAD'}", flush=True)
    if med >= 2:
        return
    feat2pt = {(int(ti[k]), int(tf[k])): int(tp[k]) for k in range(len(ti))}

    for tmin in triangle_mins:
      # ---- mutual edges -> connected components (clusters) ----
      t0 = time.time()
      pairs, counts, feat_pairs, _ = mutual_knn_matches(
        corpus, image_starts, k=K, triangle_min=tmin,
        num_trees=20, max_leaf_checks=1000, seed=1)
      pairs, counts, feat_pairs = np.asarray(pairs), np.asarray(counts), np.asarray(feat_pairs)
      imgA = np.repeat(pairs[:, 0], counts); imgB = np.repeat(pairs[:, 1], counts)
      ga = image_starts[imgA] + feat_pairs[:, 0]
      gb = image_starts[imgB] + feat_pairs[:, 1]
      nodes = np.unique(np.concatenate([ga, gb]))
      idx = {int(g): i for i, g in enumerate(nodes)}
      ia = np.fromiter((idx[int(g)] for g in ga), dtype=np.int64, count=len(ga))
      ib = np.fromiter((idx[int(g)] for g in gb), dtype=np.int64, count=len(gb))
      U = len(nodes)
      graph = coo_matrix((np.ones(len(ia)), (ia, ib)), shape=(U, U))
      n_clusters, labels = connected_components(graph, directed=False)
      t_cluster = time.time() - t0

      node_img = (np.searchsorted(image_starts, nodes, side="right") - 1).astype(np.int64)
      node_feat = (nodes - image_starts[node_img]).astype(np.int64)
      order = np.argsort(labels, kind="stable")
      m_img, m_feat, m_lab = node_img[order], node_feat[order], labels[order]
      sizes = np.bincount(labels, minlength=n_clusters)
      offsets = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)

      sz = sizes[sizes >= 2]
      print(f"\n-- triangle_min={tmin} --  edges={len(ga):,}  nodes={U:,}  "
            f"clusters={n_clusters:,}  (build {t_cluster:.2f}s)")
      print(f"   cluster sizes: median {int(np.median(sz))}  p90 {int(np.percentile(sz, 90))}  "
            f"max {int(sizes.max())}  (size-2: {int((sizes == 2).sum()):,})")

      # pass 1: triangulate all, trim outlier members
      t0 = time.time()
      _, _, err1, dep1, _ = triangulate_and_reproject(
          m_img, m_feat, offsets, rays_world, centers, kp, R, t, cams, N)
      keep1 = (err1 < REPROJ_PX) & (dep1 > 0)

      # pass 2: re-triangulate from survivors, drop clusters with <2
      s_img, s_feat, s_lab = m_img[keep1], m_feat[keep1], m_lab[keep1]
      surv_sz = np.bincount(s_lab, minlength=n_clusters)
      keep_member = (surv_sz >= 2)[s_lab]
      s_img, s_feat, s_lab = s_img[keep_member], s_feat[keep_member], s_lab[keep_member]
      o2 = np.argsort(s_lab, kind="stable")
      s_img, s_feat, s_lab = s_img[o2], s_feat[o2], s_lab[o2]
      _, comp = np.unique(s_lab, return_inverse=True)
      csz = np.bincount(comp)
      coff = np.concatenate([[0], np.cumsum(csz)]).astype(np.int64)
      _, _, err2, dep2, spread = triangulate_and_reproject(
          s_img, s_feat, coff, rays_world, centers, kp, R, t, cams, N)
      t_test = time.time() - t0

      ok_mem = (err2 < REPROJ_PX) & (dep2 > 0)
      all_ok = np.minimum.reduceat(ok_mem.astype(np.int8), coff[:-1].astype(np.intp)) == 1
      accept = all_ok & (spread >= MIN_ANGLE_DEG) & (csz >= 2)
      n_acc = int(accept.sum())

      # reconcile accepted clusters with existing tracks
      cls = Counter(); new_obs = 0; dup_img = 0
      for c in np.nonzero(accept)[0]:
          sl = slice(coff[c], coff[c + 1])
          imgs = s_img[sl]; feats = s_feat[sl]
          if len(set(imgs.tolist())) != len(imgs):
              dup_img += 1
          pids = {feat2pt[(int(i), int(f))] for i, f in zip(imgs, feats) if (int(i), int(f)) in feat2pt}
          linked = sum(1 for i, f in zip(imgs, feats) if (int(i), int(f)) in feat2pt)
          if len(pids) == 0:
              cls["new"] += 1; new_obs += len(imgs)
          elif len(pids) == 1:
              cls["extend"] += 1; new_obs += len(imgs) - linked
          else:
              cls["merge"] += 1
      avg = np.mean([coff[c + 1] - coff[c] for c in np.nonzero(accept)[0]]) if n_acc else 0
      big = np.nonzero(sizes >= 10)[0]
      bacc = np.isin(big, np.nonzero(accept)[0]).sum() if len(big) else 0

      print(f"   geom test {t_test:.2f}s (no RANSAC)  accepted clusters="
            f"{n_acc:,}/{int((sizes >= 2).sum()):,}  avg {avg:.1f} views  "
            f"dup-image {dup_img:,}")
      print(f"   candidate TRACKS: new={cls['new']:,}  extend={cls['extend']:,}  "
            f"merge={cls['merge']:,}  ~{new_obs:,} new obs   |  "
            f"clusters>=10: {len(big):,}, accepted {int(bacc):,}")


if __name__ == "__main__":
    base = Path("/tmp")
    for label, ws in [("seattle_backyard (WIN)", base / "seattle_backyard_ab_ws"),
                      ("seoul_bull (small)", base / "seoul_bull_sculpture_ab_ws")]:
        recon = ws / "cluster_recon.sfmr"
        if not recon.exists():
            print(f"skip {label}: no {recon}", file=sys.stderr); continue
        run(label, recon, next((ws / "features").glob("sift-*")),
            triangle_mins=(0, 1, 2, 4))
