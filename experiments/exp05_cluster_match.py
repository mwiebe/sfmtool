# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 05 — cluster-based pairwise matching.

Turns the descriptor-clustering pipeline (exp02–04) into a real matcher whose
output is a ``.matches`` file consumable by ``sfm solve -i``:

    1. all descriptors into one approximate-NN index (the in-tree
       ``sfmtool.KdForest`` randomized kd-tree forest; ``--exact`` for the
       brute-force oracle),
    2. data-derived bounded-radius graph -> connected components = candidate
       tracks (the auto-threshold from exp04, optionally with the
       isolated-point prefilter),
    3. for each cluster, emit all in-cluster *cross-image* pairs, bucketed by
       (image_i, image_j) -> [(feature_i, feature_j), ...],
    4. write a ``.matches`` file with ``has_two_view_geometries=False`` so
       ``sfm solve -i`` runs COLMAP's own two-view verification on the pairs.

Tuned for **recall**, not purity: we *want* a looser T than the exp03 purity
sweet spot because COLMAP's RANSAC verifier in step 4 is the precision filter.

Usage:
    pixi run -e experiments python experiments/exp05_cluster_match.py \
        ../seoul_bull_ws/sfmr/*.sfmr \
        --out ../seoul_bull_ws/matches/cluster.matches
"""

from __future__ import annotations

import argparse
import glob
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from exp03_radius_clusters import K, knn_all
from seed_cluster import seed_claim_clusters
from exp04_auto_threshold import gmm_crossover, otsu_threshold
from sfm_descriptors import load_descriptor_bank
from sfmtool._sfmtool import read_sift_metadata, write_matches


def derive_threshold(d1: np.ndarray, method: str) -> float:
    d1f = d1[(d1 > 1.0) & (d1 <= np.percentile(d1, 99.5))]
    if method == "otsu":
        return otsu_threshold(d1f)
    if method == "gmm":
        return gmm_crossover(d1f)
    return float(method)


def cliff_threshold(dst: np.ndarray, pct: float = 50.0) -> float:
    """Per-point radius-estimate view (spec §2b): T = the ``pct`` percentile of
    the first-excluded distance d_(r+1) — the neighbour just past each point's
    cliff (the largest relative jump in its sorted neighbours). Reads the shared
    k-NN distance table, not a separate query. Default pct=50 (the median)."""
    dn = dst[:, 1:8]
    g = (dn[:, 1:] / np.maximum(dn[:, :-1], 1e-6)).argmax(1)
    outer = dst[np.arange(len(dst)), 2 + g]
    return float(np.percentile(outer, pct))


def build_pairs(bank, idx, dst, T, active, max_pairs_per_cluster):
    """Cluster the descriptors and emit cross-image pair lists per image pair.

    Returns ``image_pair -> [(feat_i, feat_j, distance), ...]`` plus per-cluster
    stats (sizes, max-image-coverage), and the cluster labels (-1 for inactive
    or singletons).
    """
    n = bank.n
    img = bank.image_label
    i_rep = np.repeat(np.arange(n), K)
    j = idx.ravel().astype(np.int64)
    dd = dst.ravel()
    keep = (dd <= T) & (i_rep != j) & active[i_rep] & active[j] & (img[i_rep] != img[j])
    rows = i_rep[keep]
    cols = j[keep]
    g = csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n, n))
    _, comp = connected_components(g, directed=False)
    sizes = np.bincount(comp, minlength=n)
    # mark inactive/singleton descriptors as -1 so we don't emit pairs from them
    labels = np.where(active & (sizes[comp] >= 2), comp, -1)

    # Group members by component (only multi-component active ones).
    members_by_comp: dict[int, np.ndarray] = {}
    valid = np.flatnonzero(labels >= 0)
    if valid.size:
        order = valid[np.argsort(labels[valid], kind="stable")]
        clab = labels[order]  # sorted component labels along `order`
        bounds = np.r_[0, np.flatnonzero(np.diff(clab)) + 1, len(order)]
        for s, e in zip(bounds[:-1], bounds[1:]):
            members_by_comp[int(clab[s])] = order[s:e]

    # A real track can have at most one feature per image, so a cluster bigger
    # than the image count cannot be a single physical point — it's a mega-
    # cluster bridging repeated structure.  Drop it; let the rest of the
    # candidate clusters supply the matches.  This is a knob-free upper bound.
    n_images = int(img.max() + 1)
    dropped_mega = 0
    # Build the (image_pair) -> [(feat_i, feat_j, dist)] bucket map.
    pair_map: dict[tuple[int, int], list[tuple[int, int, float]]] = {}
    for c, members in members_by_comp.items():
        if len(members) > n_images:
            dropped_mega += 1
            continue
        if len(members) > max_pairs_per_cluster:
            members = members[:max_pairs_per_cluster]
        # all-pairs cross-image
        imgs = img[members]
        order = np.argsort(imgs, kind="stable")
        members = members[order]
        imgs = imgs[order]
        for a in range(len(members)):
            for b in range(a + 1, len(members)):
                ia, ib = int(imgs[a]), int(imgs[b])
                if ia == ib:
                    continue
                ma, mb = int(members[a]), int(members[b])
                # descriptor L2 distance between the two members
                d = float(
                    np.linalg.norm(
                        bank.descriptors[ma].astype(np.float32)
                        - bank.descriptors[mb].astype(np.float32)
                    )
                )
                fa, fb = int(bank.feature_label[ma]), int(bank.feature_label[mb])
                key = (ia, ib)  # ordered: ia < ib because we sorted by img
                pair_map.setdefault(key, []).append((fa, fb, d))
    if dropped_mega:
        print(f"dropped {dropped_mega} mega-cluster(s) (size > {n_images} images)")
    return pair_map, members_by_comp, labels


def build_pairs_neighbors(bank, idx, dst, T, active):
    """No transitive merging: each descriptor's own within-radius cross-image
    neighbours (from its 16-NN) are the match candidates.  Avoids the
    connected-component chaining / mega-clusters entirely and densifies pairs
    (every descriptor contributes, not one match per cluster).
    """
    n = bank.n
    img = bank.image_label
    feat = bank.feature_label
    i_rep = np.repeat(np.arange(n), K)
    j = idx.ravel().astype(np.int64)
    dd = dst.ravel()
    keep = (dd <= T) & (i_rep != j) & active[i_rep] & active[j] & (img[i_rep] != img[j])
    src = i_rep[keep]
    dstn = j[keep]
    dval = dd[keep]
    # Order each edge by image index so (i,j) and (j,i) collapse to one bucket.
    ai = img[src]
    aj = img[dstn]
    lo_is_src = ai < aj
    img_lo = np.where(lo_is_src, ai, aj)
    img_hi = np.where(lo_is_src, aj, ai)
    feat_lo = np.where(lo_is_src, feat[src], feat[dstn])
    feat_hi = np.where(lo_is_src, feat[dstn], feat[src])

    pair_map: dict[tuple[int, int], list[tuple[int, int, float]]] = {}
    for il, ih, fl, fh, d in zip(img_lo, img_hi, feat_lo, feat_hi, dval):
        pair_map.setdefault((int(il), int(ih)), []).append((int(fl), int(fh), float(d)))
    return pair_map


def resolve_one_per_image(
    pair_map: dict[tuple[int, int], list[tuple[int, int, float]]],
) -> dict[tuple[int, int], list[tuple[int, int, float]]]:
    """For each image pair, enforce one-to-one feature matching by greedy
    smallest-distance assignment.  This is a cheap, local one-per-image
    resolution; the verifier downstream will further enforce geometric
    one-to-one.
    """
    out: dict[tuple[int, int], list[tuple[int, int, float]]] = {}
    for pair, items in pair_map.items():
        items.sort(key=lambda x: x[2])
        used_i: set[int] = set()
        used_j: set[int] = set()
        kept: list[tuple[int, int, float]] = []
        for fa, fb, d in items:
            if fa in used_i or fb in used_j:
                continue
            used_i.add(fa)
            used_j.add(fb)
            kept.append((fa, fb, d))
        if kept:
            out[pair] = kept
    return out


def keep_min_by_key(keys: np.ndarray, vals: np.ndarray) -> np.ndarray:
    """Boolean mask selecting the minimum-`val` element for each distinct key."""
    if len(keys) == 0:
        return np.zeros(0, dtype=bool)
    order = np.lexsort((vals, keys))  # sort by key, then val ascending
    ks = keys[order]
    first = np.empty(len(ks), dtype=bool)
    first[0] = True
    first[1:] = ks[1:] != ks[:-1]
    mask = np.zeros(len(keys), dtype=bool)
    mask[order[first]] = True
    return mask


def build_neighbor_matches_arrays(bank, idx, dst, T, active):
    """Vectorised neighbours-mode matcher: each descriptor's within-radius
    cross-image 16-NN edges, deduped one-to-one per image pair (two-pass
    mutual-best), returned as the parallel arrays write_matches wants.
    """
    n = bank.n
    img = bank.image_label.astype(np.int64)
    feat = bank.feature_label.astype(np.int64)
    i_rep = np.repeat(np.arange(n), K)
    j = idx.ravel().astype(np.int64)
    dd = dst.ravel().astype(np.float32)
    keep = (dd <= T) & (i_rep != j) & active[i_rep] & active[j] & (img[i_rep] != img[j])
    src, dstn, dval = i_rep[keep], j[keep], dd[keep]
    ai, aj = img[src], img[dstn]
    lo_is_src = ai < aj
    img_lo = np.where(lo_is_src, ai, aj)
    img_hi = np.where(lo_is_src, aj, ai)
    feat_lo = np.where(lo_is_src, feat[src], feat[dstn])
    feat_hi = np.where(lo_is_src, feat[dstn], feat[src])

    n_images = int(img.max()) + 1
    maxf = int(max(feat_lo.max(initial=0), feat_hi.max(initial=0))) + 1
    pair_key = img_lo * n_images + img_hi

    # One-to-one per pair: best target per (pair, source feat), then best source
    # per (pair, target feat).
    m1 = keep_min_by_key(pair_key * maxf + feat_lo, dval)
    pair_key, feat_lo, feat_hi, dval = (
        pair_key[m1],
        feat_lo[m1],
        feat_hi[m1],
        dval[m1],
    )
    m2 = keep_min_by_key(pair_key * maxf + feat_hi, dval)
    pair_key, feat_lo, feat_hi, dval = (
        pair_key[m2],
        feat_lo[m2],
        feat_hi[m2],
        dval[m2],
    )

    # Group by pair_key into the parallel output arrays.
    order = np.argsort(pair_key, kind="stable")
    pk = pair_key[order]
    fl = feat_lo[order].astype(np.uint32)
    fh = feat_hi[order].astype(np.uint32)
    dv = dval[order].astype(np.float32)
    if len(pk):
        bounds = np.r_[0, np.flatnonzero(np.diff(pk)) + 1]
        upk = pk[bounds]
        counts = np.diff(np.r_[bounds, len(pk)])
        image_index_pairs = np.stack([upk // n_images, upk % n_images], axis=1).astype(
            np.uint32
        )
    else:
        image_index_pairs = np.zeros((0, 2), dtype=np.uint32)
        counts = np.zeros(0, dtype=np.uint32)
    match_counts = counts.astype(np.uint32)
    feat_idx = np.stack([fl, fh], axis=1)
    return image_index_pairs, match_counts, feat_idx, dv


def build_cluster_matches_arrays(
    bank, idx, dst, T, min_size=2, refine_iters=0, forest=None
):
    """Materialize density-seeded clusters and emit all in-cluster cross-image
    pairs. Clusters are a hard partition with one feature per image, so the
    resulting matches are already one-to-one per image pair. ``min_size`` keeps
    only clusters spanning at least that many images (multiplicity filter).
    ``refine_iters>0`` mean-shift refines each cluster (needs ``forest``).
    """
    labels = seed_claim_clusters(
        idx,
        dst,
        bank.image_label,
        T,
        descriptors=bank.descriptors if refine_iters else None,
        refine_iters=refine_iters,
        forest=forest,
    )
    img = bank.image_label
    feat = bank.feature_label
    X = bank.descriptors.astype(np.float32)

    valid = np.flatnonzero(labels >= 0)
    order = valid[np.argsort(labels[valid], kind="stable")]
    clab = labels[order]
    bounds = np.r_[0, np.flatnonzero(np.diff(clab)) + 1, len(order)]

    pair_map: dict[tuple[int, int], list[tuple[int, int, float]]] = {}
    for s, e in zip(bounds[:-1], bounds[1:]):
        members = order[s:e]
        if len(members) < min_size:
            continue
        imgs = img[members]
        o = np.argsort(imgs, kind="stable")
        members, imgs = members[o], imgs[o]
        for a in range(len(members)):
            for b in range(a + 1, len(members)):
                ia, ib = int(imgs[a]), int(imgs[b])
                if ia == ib:
                    continue
                ma, mb = int(members[a]), int(members[b])
                d = float(np.linalg.norm(X[ma] - X[mb]))
                pair_map.setdefault((ia, ib), []).append(
                    (int(feat[ma]), int(feat[mb]), d)
                )
    n_clusters = int((labels >= 0).any()) and int(labels.max()) + 1
    print(f"materialized {n_clusters} clusters")
    return pairmap_to_arrays(pair_map)


def pairmap_to_arrays(pair_map):
    """Convert the components-mode dict to the parallel output arrays."""
    sorted_pairs = sorted(pair_map.keys())
    image_index_pairs = (
        np.array(sorted_pairs, dtype=np.uint32)
        if sorted_pairs
        else np.zeros((0, 2), dtype=np.uint32)
    )
    match_counts = np.array([len(pair_map[p]) for p in sorted_pairs], dtype=np.uint32)
    total = int(match_counts.sum())
    feat_idx = np.empty((total, 2), dtype=np.uint32)
    dists = np.empty(total, dtype=np.float32)
    off = 0
    for p in sorted_pairs:
        rows = pair_map[p]
        m = len(rows)
        feat_idx[off : off + m] = np.array(
            [(a, b) for a, b, _ in rows], dtype=np.uint32
        )
        dists[off : off + m] = np.array([d for _, _, d in rows], dtype=np.float32)
        off += m
    return image_index_pairs, match_counts, feat_idx, dists


def assemble_matches_dict(
    bank,
    image_index_pairs: np.ndarray,
    match_counts: np.ndarray,
    feat_idx: np.ndarray,
    dists: np.ndarray,
    image_names: list[str],
    sift_paths: list[Path],
    workspace_dir: Path,
    output_path: Path,
) -> dict:
    """Build the matches_data dict expected by ``write_matches``.

    Takes the already-assembled parallel arrays:
      image_index_pairs (P,2) uint32, match_counts (P,) uint32,
      feat_idx (M,2) uint32, dists (M,) float32.
    """
    total = int(match_counts.sum())

    # Per-image SIFT content & feature-tool hashes (workspace integrity check).
    feature_tool_hashes = np.zeros((len(image_names), 16), dtype=np.uint8)
    sift_content_hashes = np.zeros((len(image_names), 16), dtype=np.uint8)
    feature_counts = np.zeros(len(image_names), dtype=np.uint32)
    for i, sp in enumerate(sift_paths):
        meta = read_sift_metadata(str(sp))
        ch = meta["content_hash"]
        feature_tool_hashes[i] = np.frombuffer(
            bytes.fromhex(ch["feature_tool_xxh128"]), dtype=np.uint8
        )
        sift_content_hashes[i] = np.frombuffer(
            bytes.fromhex(ch["content_xxh128"]), dtype=np.uint8
        )
        feature_counts[i] = int(meta["metadata"]["feature_count"])

    out_abs = Path(os.path.abspath(output_path))
    workspace_rel = os.path.relpath(workspace_dir, out_abs.parent).replace("\\", "/")
    # Pull the workspace's actual feature config (feature_tool, feature_type,
    # feature_options, feature_prefix_dir) so `sfm solve -i` resolves .sift
    # files at the same path the workspace stores them.
    import json

    ws_cfg = json.loads((workspace_dir / ".sfm-workspace.json").read_text())
    return {
        "metadata": {
            "version": 1,
            "matching_method": "exhaustive",
            "matching_tool": "sfmtool-cluster-poc",
            "matching_tool_version": "",
            "matching_options": {},
            "workspace": {
                "absolute_path": str(workspace_dir),
                "relative_path": workspace_rel,
                "contents": {
                    "feature_tool": ws_cfg.get("feature_tool", "colmap"),
                    "feature_type": ws_cfg.get("feature_type", "sift"),
                    "feature_options": ws_cfg.get("feature_options", {}) or {},
                    "feature_prefix_dir": ws_cfg.get("feature_prefix_dir", "") or "",
                },
            },
            "timestamp": datetime.now(timezone.utc).astimezone().isoformat(),
            "image_count": len(image_names),
            "image_pair_count": int(len(image_index_pairs)),
            "match_count": total,
            "has_two_view_geometries": False,
        },
        "content_hash": {},
        "image_names": image_names,
        "feature_tool_hashes": feature_tool_hashes,
        "sift_content_hashes": sift_content_hashes,
        "feature_counts": feature_counts,
        "image_index_pairs": image_index_pairs,
        "match_counts": match_counts,
        "match_feature_indexes": feat_idx,
        "match_descriptor_distances": dists,
        "has_two_view_geometries": False,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("sfmr")
    ap.add_argument("--out", required=True, help="output .matches path")
    ap.add_argument("--cache", default=None, help="17-NN cache .npz")
    ap.add_argument(
        "--threshold",
        default="cliff",
        help="cliff (p50 of first-excluded; default) / otsu / gmm / float",
    )
    ap.add_argument(
        "--cliff-pct",
        type=float,
        default=50.0,
        help="percentile for the cliff threshold (--threshold cliff). Default 50.",
    )
    ap.add_argument(
        "--t-scale",
        type=float,
        default=1.0,
        help="multiply derived T by this factor (recall-first; downstream RANSAC"
        " is the precision filter). Default 1.0; ~1.25 pairs well with otsu/gmm.",
    )
    ap.add_argument(
        "--prefilter",
        action="store_true",
        help="drop isolated points (d1>antimode or ratio>0.85)",
    )
    ap.add_argument("--max-pairs-per-cluster", type=int, default=64)
    ap.add_argument(
        "--mode",
        choices=["components", "neighbors", "clusters"],
        default="neighbors",
        help="components: connected-component clusters; "
        "neighbors: per-descriptor 16-NN neighbourhoods (no merging); "
        "clusters: materialized density-seeded clusters",
    )
    ap.add_argument(
        "--exact",
        action="store_true",
        help="use exact NN (oracle) instead of the KdForest index",
    )
    ap.add_argument(
        "--preset",
        default="accurate",
        help="KdForest preset: accurate / balanced / fast",
    )
    ap.add_argument(
        "--min-cluster-size",
        type=int,
        default=2,
        help="clusters mode: keep clusters spanning >= this many images",
    )
    ap.add_argument(
        "--refine",
        type=int,
        default=0,
        help="clusters mode: mean-shift refinement iterations",
    )
    args = ap.parse_args()

    path = sorted(glob.glob(args.sfmr))[0]
    bank = load_descriptor_bank(path)
    print(f"loaded {path}: n={bank.n}")

    forest = None
    if args.exact:
        # Oracle path: exact NN, cached.
        cache = Path(args.cache or f"out/{Path(bank.workspace_dir).name}_knn17.npz")
        if cache.exists():
            z = np.load(cache)
            idx, dst = z["idx"], z["dst"]
            print(f"exact NN (oracle), cached: {cache}")
        else:
            print(f"computing exact {K}-NN over {bank.n} (one-time) -> {cache}")
            idx, dst = knn_all(bank.descriptors, K)
            cache.parent.mkdir(parents=True, exist_ok=True)
            np.savez(cache, idx=idx, dst=dst)
    else:
        # Production path: the in-tree randomized kd-tree forest (sfmtool.KdForest).
        from sfmtool import KdForest

        print(f"building KdForest (preset={args.preset}) over {bank.n} descriptors")
        forest = KdForest(np.ascontiguousarray(bank.descriptors), preset=args.preset)
        idx, dst = forest.query(np.ascontiguousarray(bank.descriptors), k=K)
        idx = idx.astype(np.int64)

    d1 = dst[:, 1]
    d5 = dst[:, 5]
    if args.threshold == "cliff":
        T_base = cliff_threshold(dst, args.cliff_pct)
        method = f"cliff p{args.cliff_pct:g}"
    else:
        T_base = derive_threshold(d1, args.threshold)
        method = args.threshold
    T = T_base * args.t_scale
    print(
        f"threshold T_base={T_base:.1f}  -> T={T:.1f} "
        f"(scale {args.t_scale})  via {method}"
    )

    if args.prefilter:
        active = ~((d1 / np.maximum(d5, 1e-6) > 0.85) | (d1 > T_base))
        print(f"prefilter keeps {active.mean():.1%}")
    else:
        active = np.ones(bank.n, dtype=bool)

    if args.mode == "components":
        pair_map, members_by_comp, _ = build_pairs(
            bank, idx, dst, T, active, args.max_pairs_per_cluster
        )
        print(f"clusters with >=2 members: {len(members_by_comp)}")
        pair_map = resolve_one_per_image(pair_map)
        ii_pairs, m_counts, feat_idx, dists = pairmap_to_arrays(pair_map)
    elif args.mode == "clusters":
        if args.refine and forest is None:
            print(
                "warning: --refine needs the index; --exact has none, "
                "falling back to seed-only"
            )
        ii_pairs, m_counts, feat_idx, dists = build_cluster_matches_arrays(
            bank,
            idx,
            dst,
            T,
            min_size=args.min_cluster_size,
            refine_iters=args.refine,
            forest=forest,
        )
    else:
        print("mode=neighbors (no transitive merging, vectorised)")
        ii_pairs, m_counts, feat_idx, dists = build_neighbor_matches_arrays(
            bank, idx, dst, T, active
        )
    print(
        f"after one-per-image: {int(m_counts.sum())} matches across "
        f"{len(ii_pairs)} image pairs"
    )

    # Resolve image names -> .sift paths in the workspace.
    from sfmtool.sift.file import get_sift_path_for_image

    ws = Path(bank.workspace_dir)
    image_names = list(bank.image_names)
    sift_paths = [Path(get_sift_path_for_image(str(ws / nm))) for nm in image_names]

    out = Path(args.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    data = assemble_matches_dict(
        bank, ii_pairs, m_counts, feat_idx, dists, image_names, sift_paths, ws, out
    )
    write_matches(str(out), data)
    print(
        f"\nwrote {out}: {data['metadata']['image_pair_count']} pairs, "
        f"{data['metadata']['match_count']} matches "
        f"(has_two_view_geometries=False)"
    )
    print(
        "next: `pixi run sfm solve -i <matches>` will run COLMAP "
        "verification on these pairs."
    )


if __name__ == "__main__":
    main()
