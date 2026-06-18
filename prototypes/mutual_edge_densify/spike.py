"""Spike: mutual-kNN edges as post-reconstruction track-densification candidates.

Hypothesis (from specs/core/mutual-knn-matching.md, "Opportunity"): the raw
mutual edges are cheap to compute on top of a cluster run, and once a
reconstruction exists each edge can be validated against *known* geometry far
more cheaply than blind two-view RANSAC -- and non-overlapping pairs fall out
for free.

This is a measurement spike, not a feature. It loads a cluster-based
reconstruction, self-validates the pose math against existing tracks, computes
the raw mutual edge set, and geometrically tests every edge by triangulating
from the known camera centers and reprojecting. It reports: pass rate, the
densification yield (new tracks / extended tracks), the free gating of
non-overlapping image pairs, and the cost vs the standalone RANSAC verification.

Run: pixi run -e test python prototypes/mutual_edge_densify/spike.py
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
from sfmtool._sfmtool import SfmrReconstruction, mutual_knn_matches, read_sift

REPROJ_PX = 4.0      # max per-view reprojection error to accept a candidate
MIN_ANGLE_DEG = 1.5  # min triangulation angle (reject near-degenerate rays)
K = 20               # mutual-kNN k (matches the tuned default)


def quat_wxyz_to_R(q):  # world->cam rotation, COLMAP convention
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ])


def run(name: str, recon_path: Path, features_dir: Path):
    print(f"\n========== {name} ==========", flush=True)
    recon = SfmrReconstruction.load(str(recon_path))
    names = recon.image_names
    N = len(names)
    quats = recon.quaternions_wxyz
    trans = recon.translations
    cam_idx = recon.camera_indexes
    cameras = recon.cameras
    positions = recon.positions  # (M,3) existing 3D points
    ti = recon.track_image_indexes
    tf = recon.track_feature_indexes
    tp = recon.track_point_ids

    R = np.stack([quat_wxyz_to_R(q) for q in quats])              # (N,3,3) world->cam
    t = np.asarray(trans)                                         # (N,3)
    centers = np.einsum("nij,nj->ni", R.transpose(0, 2, 1), -t)  # (N,3) world cam centers
    cams = [cameras[int(c)] for c in cam_idx]

    # --- load sift corpus in recon image order (corpus idx == recon idx) ---
    kp, descs, starts = [], [], [0]
    for nm in names:
        d = read_sift(str(features_dir / f"{Path(nm).name}.sift"))
        kp.append(np.asarray(d["positions_xy"], dtype=np.float64))
        descs.append(np.asarray(d["descriptors"], dtype=np.uint8))
        starts.append(starts[-1] + len(descs[-1]))
    corpus = np.concatenate(descs, axis=0)
    image_starts = np.asarray(starts, dtype=np.uint32)

    # precompute world-space unit ray per keypoint per image
    rays_world = []
    for n in range(N):
        cray = np.asarray(cams[n].pixel_to_ray_batch(kp[n]))     # (Kn,3) cam frame
        rays_world.append((R[n].T @ cray.T).T)                   # (Kn,3) world

    # ---- self-check: reproject existing tracks, expect small error ----
    errs = []
    for n in range(N):
        m = ti == n
        if not m.any():
            continue
        X = positions[tp[m]]                                     # (k,3)
        xc = (R[n] @ X.T).T + t[n]
        pred = np.asarray(cams[n].ray_to_pixel_batch(xc))
        obs = kp[n][tf[m]]
        e = np.linalg.norm(pred - obs, axis=1)
        errs.append(e[np.isfinite(e)])
    errs = np.concatenate(errs)
    med = float(np.median(errs))
    print(f"pose self-check: {len(errs):,} existing obs, median reproj "
          f"{med:.3f}px, p90 {np.percentile(errs, 90):.3f}px "
          f"-> {'OK' if med < 2 else 'BAD POSE MATH'}", flush=True)
    if med >= 2:
        print("  aborting: pose/convention wrong, downstream numbers untrustworthy")
        return

    # feature -> existing point id map
    feat2pt = {(int(ti[k]), int(tf[k])): int(tp[k]) for k in range(len(ti))}

    # covisibility (overlapping image pairs share >=1 existing point)
    from collections import defaultdict
    pt_imgs = defaultdict(set)
    for k in range(len(ti)):
        pt_imgs[int(tp[k])].add(int(ti[k]))
    covis = set()
    for imgs in pt_imgs.values():
        il = sorted(imgs)
        for a in range(len(il)):
            for b in range(a + 1, len(il)):
                covis.add((il[a], il[b]))

    # ---- raw mutual edges (pre-verification) ----
    t0 = time.time()
    pairs, counts, feat_pairs, _dist = mutual_knn_matches(
        corpus, image_starts, k=K, num_trees=20, max_leaf_checks=1000, seed=1)
    t_edges = time.time() - t0
    pairs = np.asarray(pairs); counts = np.asarray(counts)
    feat_pairs = np.asarray(feat_pairs)
    M = len(feat_pairs)

    # expand per-candidate image indices from the per-pair CSR
    imgA = np.repeat(pairs[:, 0], counts)
    imgB = np.repeat(pairs[:, 1], counts)
    featA = feat_pairs[:, 0]; featB = feat_pairs[:, 1]

    # ---- geometric test: triangulate from known centers, reproject ----
    t0 = time.time()
    dirsA = np.stack([rays_world[i][f] for i, f in zip(imgA, featA)])
    dirsB = np.stack([rays_world[i][f] for i, f in zip(imgB, featB)])
    dirs = np.empty((2 * M, 3)); dirs[0::2] = dirsA; dirs[1::2] = dirsB
    cen = np.empty((2 * M, 3)); cen[0::2] = centers[imgA]; cen[1::2] = centers[imgB]
    offsets = np.arange(0, 2 * M + 1, 2, dtype=np.int64)
    from sfmtool._sfmtool import triangulate_batch
    tri = triangulate_batch(dirs, cen, offsets)
    X = np.asarray(tri["points"]); in_front = np.asarray(tri["in_front_of_all_cameras"])

    # triangulation angle between the two rays
    cosang = np.clip(np.sum(dirsA * dirsB, axis=1), -1, 1)
    angle_ok = np.degrees(np.arccos(cosang)) >= MIN_ANGLE_DEG

    # reprojection error in both views (group by image)
    errA = np.full(M, np.inf); errB = np.full(M, np.inf)
    pxA = np.stack([kp[i][f] for i, f in zip(imgA, featA)])
    pxB = np.stack([kp[i][f] for i, f in zip(imgB, featB)])
    for n in range(N):
        ma = imgA == n
        if ma.any():
            xc = (R[n] @ X[ma].T).T + t[n]
            pr = np.asarray(cams[n].ray_to_pixel_batch(xc))
            errA[ma] = np.linalg.norm(pr - pxA[ma], axis=1)
        mb = imgB == n
        if mb.any():
            xc = (R[n] @ X[mb].T).T + t[n]
            pr = np.asarray(cams[n].ray_to_pixel_batch(xc))
            errB[mb] = np.linalg.norm(pr - pxB[mb], axis=1)
    reproj_ok = np.fmax(errA, errB) < REPROJ_PX
    accept = in_front & angle_ok & reproj_ok
    t_test = time.time() - t0

    # ---- accounting ----
    n_pairs = len(pairs)
    n_covis_pairs = sum(1 for p in pairs if (int(p[0]), int(p[1])) in covis)
    pair_is_covis = np.array([(int(imgA[i]), int(imgB[i])) in covis for i in range(M)])

    # classify accepted candidates by track membership
    cls = {"both_same": 0, "merge": 0, "extend": 0, "new": 0}
    for i in np.nonzero(accept)[0]:
        a = feat2pt.get((int(imgA[i]), int(featA[i])))
        b = feat2pt.get((int(imgB[i]), int(featB[i])))
        if a is not None and b is not None:
            cls["both_same" if a == b else "merge"] += 1
        elif a is not None or b is not None:
            cls["extend"] += 1
        else:
            cls["new"] += 1

    print(f"images={N}  candidates={M:,} on {n_pairs} pairs "
          f"({n_covis_pairs} overlapping / {n_pairs - n_covis_pairs} not)")
    print(f"edge compute={t_edges:.2f}s   geometric test={t_test:.2f}s "
          f"(no RANSAC) for {M:,} candidates")
    print(f"accepted={accept.sum():,} ({100 * accept.mean():.1f}%)   "
          f"in_front={100 * in_front.mean():.0f}%  angle_ok={100 * angle_ok.mean():.0f}%  "
          f"reproj_ok={100 * reproj_ok.mean():.0f}%")
    # gating: pass rate on overlapping vs non-overlapping pairs
    for label, mask in (("overlapping", pair_is_covis), ("non-overlap", ~pair_is_covis)):
        if mask.any():
            print(f"  {label:11s}: {mask.sum():>7,} cand  accept "
                  f"{100 * accept[mask].mean():4.1f}%")
    print(f"densification yield (accepted): new-track={cls['new']:,}  "
          f"extend-track={cls['extend']:,}  merge={cls['merge']:,}  "
          f"already-linked={cls['both_same']:,}")


if __name__ == "__main__":
    base = Path("/tmp")
    datasets = [
        ("seattle_backyard (WIN case)", base / "seattle_backyard_ab_ws"),
        ("seoul_bull (small)", base / "seoul_bull_sculpture_ab_ws"),
    ]
    for label, ws in datasets:
        recon = ws / "cluster_recon.sfmr"
        feats = next((ws / "features").glob("sift-*"))
        if not recon.exists():
            print(f"skip {label}: no {recon}", file=sys.stderr)
            continue
        run(label, recon, feats)
