# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment: the seed's opening moves, written against the core APIs alone.

A deliberately small re-expression of rung 1's first steps -- base admission,
coarse cut, capture intrinsics, covisibility seed groups, and one ray-space
probe per group -- as an evaluation of the sfmtool core surface: does each
idea come out as one named call, and where does numpy glue creep in?

Read-only against the workspace: opens the given cluster ``.matches`` file and
prints a transcript.  When ``SFMTOOL_OPENING_METRICS`` names a path, the
per-proposal metric rows are additionally written there as JSON (the input for
the classification-utility mining); nothing else is ever written.

    pixi run -e dev python scripts/exp_seed_opening.py <clusters.matches>
"""

import json
import os
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

from sfmtool._sfmtool.analysis import coarsest_cluster_ids, triangulate_batch
from sfmtool._sfmtool.geometry import estimate_essential_rays, estimate_intrinsics
from sfmtool._sfmtool.io import MatchesFile
from sfmtool._sfmtool.matching import ClusterCovisibility

N_COARSE = 3000  # the fleet's seedability floor (see exp_fast_seed.py)
GROUP_SIZE = 5  # images per covisibility seed group
N_PROBES = 8  # probe outcomes the drive loop pulls proposals until it holds
TOL_PX = 3.0  # keypoint localization tolerance behind the ray consensus bound
SEED = 0


def member_arrays(sel):
    """A selection's member-parallel view: (cluster, image, position) rows.

    API-fit note: the file stores `cluster_starts`, and the intrinsics vote
    now takes the selection itself; the one consumer left for the expanded
    CSR index is this script's own pair join (step 5), where the
    member-parallel cluster column is the working shape.  That consumer has
    shrunk to one join per probed group, since a proposal names its seed
    pair -- the expansion survives for the join's shape, not its volume.
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
    Called once per probed group, on the pair the proposal already named.
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
    """Ray-space two-view init: the pair's arbitration metrics, as a dict.

    The ray-native reading throughout: depth positive ALONG THE RAY in both
    cameras (`z > 0` is wrong past 90 degrees), midpoint triangulation.
    Returns None when no consensus essential exists or cheiral support is
    under 12 rays; otherwise n_corr / n_inlier / cheiral counts, median and
    p90 triangulation parallax, the relative rotation magnitude, and the
    estimator's essentialness -- the candidate features for classifying the
    group's motion regime (a pan reads large rot_deg at near-zero parallax).
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
        if best is None or good.sum() > best[1].sum():
            best = (rot, good, x, -rot.T @ t)
    rot, good, x, c2 = best
    if good.sum() < 12:
        return None
    xa = x[good]
    v1 = xa / np.linalg.norm(xa, axis=1, keepdims=True)
    v2 = xa - c2
    v2 /= np.linalg.norm(v2, axis=1, keepdims=True)
    par = np.degrees(np.arccos(np.clip((v1 * v2).sum(1), -1, 1)))
    return {
        "n_corr": int(len(x1)),
        "n_inlier": int(inl.sum()),
        "cheiral": int(good.sum()),
        "par_med": float(np.median(par)),
        "par_p90": float(np.percentile(par, 90)),
        "rot_deg": float(
            np.degrees(np.arccos(np.clip((np.trace(rot) - 1) / 2, -1, 1)))
        ),
        "essentialness": float(est["essentialness"]),
    }


def main():
    path = Path(sys.argv[1])

    # 1. Base admission: the loader's predicate as one native derivation.
    matches = MatchesFile(path)
    sel = matches.select_clusters(min_span=2)
    print(
        f"{path.name}: {sel.image_count} images "
        f"{sel.image_width}x{sel.image_height}, "
        f"{sel.cluster_count} clusters, {sel.member_count} members"
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
    keep = coarsest_cluster_ids(sel, N_COARSE)
    if len(keep) < sel.cluster_count:
        total = sel.cluster_count
        sel = sel.select_clusters(restrict_cluster_ids=keep)
        print(f"coarse admission: kept {len(keep)}/{total} coarsest clusters")

    # 4. Covisibility seed groups, off the COARSE working set: the proposals
    #    are read from the same graph the probe judges them on, so a group's
    #    reported edge weight is the join the probe will actually get.  (The
    #    referee principle of step 2 -- the intrinsics vote measured on the
    #    full admission -- is a different role and stays there.)
    covis = ClusterCovisibility.from_matches(sel)

    # 5. The drive loop pulls proposals while it still wants probe outcomes,
    #    and probes TWO pairs per proposal chosen for opposite purposes: the
    #    seed pair (maximum shared clusters -- best-conditioned essential) and
    #    the widest pair (maximum mean keypoint displacement -- the honest
    #    parallax ceiling; the max-shared pair of a video group is its most
    #    adjacent, so it reads parallax near the group's minimum).  A starved
    #    proposal is skipped joinlessly off seed_shared and costs no budget.
    #    Each proposal's metric row goes to SFMTOOL_OPENING_METRICS when set.
    obs_c, obs_i, uv = member_arrays(sel)
    probed = starved = 0
    rows = []
    for prop in covis.seed_image_groups(GROUP_SIZE):
        if probed >= N_PROBES:
            break
        # `images` is a uint32 numpy array; the transcript prints the group as
        # a plain list, so it is converted for the print rather than formatted
        # through numpy's own repr.
        group = prop.images.tolist()
        (a, b), shared = prop.seed_pair, prop.seed_shared
        pair_shared = prop.pair_shared
        disp = prop.pair_displacement
        pairs = list(combinations(range(len(group)), 2))  # condensed order
        row = {
            "images": group,
            "seed_pair": [int(a), int(b)],
            "seed_shared": int(shared),
            "pair_shared_min": int(pair_shared.min()),
            "pair_shared_med": float(np.median(pair_shared)),
            "disp_med": float(np.median(disp)) if disp is not None else None,
            "disp_max": float(disp.max()) if disp is not None else None,
        }
        rows.append(row)
        if shared < 20:
            starved += 1
            row["starved"] = True
            print(f"  group {group}: starved ({shared} shared clusters at best)")
            continue
        probed += 1
        x1, x2 = pair_correspondences(obs_c, obs_i, uv, a, b)
        row["seed_probe"] = seed_probe = probe_pair(cam, x1, x2, f)
        if seed_probe is None:
            print(f"  group {group}: seed ({a}, {b}) no consensus")
        else:
            print(
                f"  group {group}: seed ({a}, {b}) {seed_probe['cheiral']} cheiral, "
                f"parallax {seed_probe['par_med']:.2f}|{seed_probe['par_p90']:.2f} deg, "
                f"rot {seed_probe['rot_deg']:.2f} deg, "
                f"essentialness {seed_probe['essentialness']:.4f}"
            )
        # The widest pair, when it is a different pair with enough shared
        # clusters to estimate on at all.
        if disp is not None:
            wa, wb = (group[i] for i in pairs[int(np.argmax(disp))])
            w_shared = int(pair_shared[int(np.argmax(disp))])
            row["wide_pair"] = [int(wa), int(wb)]
            row["wide_shared"] = w_shared
            row["wide_disp_px"] = float(disp.max())
            if (wa, wb) != (a, b) and w_shared >= 8:
                x1, x2 = pair_correspondences(obs_c, obs_i, uv, wa, wb)
                row["wide_probe"] = wide = probe_pair(cam, x1, x2, f)
                if wide is None:
                    print(f"    wide ({wa}, {wb}) shared {w_shared}: no consensus")
                else:
                    print(
                        f"    wide ({wa}, {wb}) shared {w_shared}, "
                        f"disp {row['wide_disp_px']:.1f}px: {wide['cheiral']} cheiral, "
                        f"parallax {wide['par_med']:.2f}|{wide['par_p90']:.2f} deg, "
                        f"rot {wide['rot_deg']:.2f} deg"
                    )
    print(f"covisibility: {probed} proposals probed, {starved} starved skipped")
    if out := os.environ.get("SFMTOOL_OPENING_METRICS"):
        Path(out).write_text(
            json.dumps({"file": path.name, "proposals": rows}, indent=1)
        )
        print(f"metrics: {len(rows)} proposal rows -> {out}")


if __name__ == "__main__":
    main()
