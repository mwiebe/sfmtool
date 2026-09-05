# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment: the seed's opening moves, written against the core APIs alone.

A deliberately small re-expression of rung 1's first steps -- base admission,
coarse cut, capture intrinsics, covisibility seed groups, and one ray-space
probe per group -- as an evaluation of the sfmtool core surface: does each
idea come out as one named call, and where does numpy glue creep in?

Read-only: opens the given cluster ``.matches`` file and prints a transcript;
nothing is written.

    pixi run -e dev python scripts/exp_seed_opening.py <clusters.matches>
"""

import sys
from itertools import combinations, islice
from pathlib import Path

import numpy as np

from sfmtool._sfmtool.analysis import coarsest_clusters, triangulate_batch
from sfmtool._sfmtool.geometry import estimate_essential_rays, estimate_intrinsics
from sfmtool._sfmtool.io import MatchesFile
from sfmtool._sfmtool.matching import ClusterCovisibility

N_COARSE = 3000  # the fleet's seedability floor (see exp_fast_seed.py)
GROUP_SIZE = 5  # images per covisibility seed group
N_GROUPS = 8  # proposals probed per run
TOL_PX = 3.0  # keypoint localization tolerance behind the ray consensus bound
SEED = 0


def member_arrays(sel):
    """A selection's member-parallel view: (cluster, image, position) rows.

    API-fit note: the file stores `cluster_starts`, and the intrinsics vote
    now takes the selection itself; the one consumer left for the expanded
    CSR index is this script's own pair joins (step 5), where the
    member-parallel cluster column is the working shape.
    """
    starts = np.asarray(sel.cluster_starts, dtype=np.int64)
    obs_c = np.repeat(np.arange(len(starts) - 1), np.diff(starts))
    obs_i = np.asarray(sel.member_images, dtype=np.int64)
    uv = np.asarray(sel.member_positions(), dtype=np.float64)
    return obs_c, obs_i, uv


def pair_correspondences(obs_c, obs_i, uv, img_a, img_b):
    """Positions of the clusters two images share, one row per cluster.

    API-fit note: numpy glue.  A selection keeps at most one reference/kept
    member per (cluster, image), so a cluster->row map per image joins them.
    """
    n_cl = obs_c.max() + 1
    row = {img: np.full(n_cl, -1, dtype=np.int64) for img in (img_a, img_b)}
    for img in (img_a, img_b):
        mask = obs_i == img
        row[img][obs_c[mask]] = np.nonzero(mask)[0]
    shared = (row[img_a] >= 0) & (row[img_b] >= 0)
    return uv[row[img_a][shared]], uv[row[img_b][shared]]


def decompose_essential(e):
    """The four (R, t) readings of an essential matrix.

    API-fit note: textbook glue (U W V^T with sign fixes); the kernel returns
    `e_matrix` and leaves the decomposition to the caller.
    """
    u, _, vt = np.linalg.svd(e)
    if np.linalg.det(u) < 0:
        u = -u
    if np.linalg.det(vt) < 0:
        vt = -vt
    w = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    for rot in (u @ w @ vt, u @ w.T @ vt):
        for t in (u[:, 2], -u[:, 2]):
            yield rot, t


def probe_pair(cam, x1, x2, f):
    """Ray-space two-view init: cheiral support and median parallax.

    The ray-native reading throughout: depth positive ALONG THE RAY in both
    cameras (`z > 0` is wrong past 90 degrees), midpoint triangulation.
    API-fit note: one estimator call plus ~30 lines of decomposition,
    cheirality and parallax -- the block is a core-primitive candidate.
    """
    d1 = np.ascontiguousarray(cam.pixel_to_ray_batch(np.ascontiguousarray(x1)))
    d2 = np.ascontiguousarray(cam.pixel_to_ray_batch(np.ascontiguousarray(x2)))
    est = estimate_essential_rays(d1, d2, max_angle_rad=TOL_PX / f, seed=SEED)
    if est is None:
        return None
    inl = np.asarray(est["inliers"], dtype=bool)
    a1, a2 = np.ascontiguousarray(d1[inl]), np.ascontiguousarray(d2[inl])
    best = None
    for rot, t in decompose_essential(np.asarray(est["e_matrix"])):
        n = len(a1)
        dirs = np.empty((2 * n, 3))
        dirs[0::2], dirs[1::2] = a1, a2 @ rot  # camera-2 rays in frame 1
        centers = np.zeros((2 * n, 3))
        centers[1::2] = -rot.T @ t
        offsets = np.arange(0, 2 * n + 1, 2, dtype=np.int64)
        x = np.asarray(triangulate_batch(dirs, centers, offsets)["points"])
        good = np.isfinite(x).all(axis=1)
        good &= np.einsum("ij,ij->i", x, a1) > 0.0
        good &= np.einsum("ij,ij->i", x @ rot.T + t, a2) > 0.0
        if best is None or good.sum() > best[0].sum():
            best = (good, x, -rot.T @ t)
    good, x, c2 = best
    if good.sum() < 12:
        return None
    xa = x[good]
    v1 = xa / np.linalg.norm(xa, axis=1, keepdims=True)
    v2 = xa - c2
    v2 /= np.linalg.norm(v2, axis=1, keepdims=True)
    parallax = np.degrees(np.median(np.arccos(np.clip((v1 * v2).sum(1), -1, 1))))
    return int(good.sum()), float(parallax), float(est["essentialness"])


def main():
    path = Path(sys.argv[1])

    # 1. Base admission: the loader's predicate as one native derivation.
    matches = MatchesFile(path)
    sel = matches.select_clusters(min_span=2)
    n_cl, n_img = len(sel.cluster_starts) - 1, len(sel.image_names)
    (w, h) = sel.image_dims[0]
    print(
        f"{path.name}: {n_img} images {w}x{h}, {n_cl} clusters, "
        f"{len(sel.member_images)} members"
    )

    # 2. Capture intrinsics, on the FULL admission (the referee must not
    #    shrink with the solve's working set): focal AND camera model, one call
    #    over the selection itself -- the file already states the observations
    #    in the layout the vote takes.
    est = estimate_intrinsics(sel, seed=SEED, columns="auto")
    cam = est["camera"]
    if cam is None:
        sys.exit(f"{path.name}: no consensus focal, so no camera to probe with")
    f = cam.focal_lengths[0]
    chart = "equidistant" if cam.model == "EQUIDISTANT_FISHEYE" else "pinhole"
    vote = est["screening_vote"] or est["vote"]
    print(
        f"intrinsics: {chart} f={f:.1f} px (vote {vote['focal_px']:.1f} over "
        f"{vote['n_pool']} votes; verdict {est['camera_model']}, "
        f"escalation {est['escalation'] or 'none'})"
    )

    # 3. Coarse cut: the alias-free working set, as a derived selection.
    keep = coarsest_clusters(sel, N_COARSE)
    if len(keep) < n_cl:
        sel = sel.select_clusters(
            min_span=2, restrict_cluster_ids=[int(c) for c in keep]
        )
        print(f"coarse admission: kept {len(keep)}/{n_cl} coarsest clusters")

    # 4. Covisibility seed groups, off the capture-level graph.
    covis = ClusterCovisibility.from_matches(matches)
    groups = list(islice(covis.seed_groups(GROUP_SIZE, 8), N_GROUPS))
    print(f"covisibility: {len(groups)} seed groups proposed")

    # 5. One probe per group: the geometry arbitrates the proposals.  The
    #    pair joins are the first (and only) step that works member-parallel.
    obs_c, obs_i, uv = member_arrays(sel)
    for group in groups:
        pairs = [
            (len(pair_correspondences(obs_c, obs_i, uv, a, b)[0]), a, b)
            for a, b in combinations(group, 2)
        ]
        shared, a, b = max(pairs)
        if shared < 20:
            print(f"  group {group}: starved ({shared} shared clusters at best)")
            continue
        x1, x2 = pair_correspondences(obs_c, obs_i, uv, a, b)
        probe = probe_pair(cam, x1, x2, f)
        if probe is None:
            print(f"  group {group}: pair ({a}, {b}) no consensus")
            continue
        cheiral, parallax, essentialness = probe
        print(
            f"  group {group}: pair ({a}, {b}) {cheiral} cheiral, "
            f"parallax {parallax:.2f} deg, essentialness {essentialness:.4f}"
        )


if __name__ == "__main__":
    main()
