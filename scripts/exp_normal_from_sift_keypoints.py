#!/usr/bin/env python3
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment: initial normal refinement — reprojected centers vs. SIFT-solve
keypoints.

Production's step-1 `refine_normals` positions each track view's patch at the
*reprojection* of the triangulated 3D point. But the SIFT solve already has the
*actual detected keypoint* for every track observation; the reprojection differs
from it by the triangulation/pose residual. This puts the two head-to-head: run
the same coarse-to-fine normal search over the track twice —

  * proj: each view's patch at the reprojected center (= production)
  * sift: each view's patch shifted in-plane to land on its SIFT keypoint

— and compare the normal, the consensus Phi, and leave-one-out ZNCC. The gain
scales with the solve's reprojection error: nil on tight solves, measurable on
harder ones (see reports/exp/2026-06-21-mvs-normal-refinement.md).

Usage::

    pixi run python scripts/exp_normal_from_sift_keypoints.py RECON.sfmr [N]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _normal_strip_lib import (  # noqa: E402
    ANG_RANGE,
    EXTENT_FACTOR,
    PATCH,
    gauss_window,
)
from exp_reference_refine import loo_zncc  # noqa: E402

# Reuse the normal-search machinery from the companion keypoints experiment.
from exp_normal_from_congealed_keypoints import (  # noqa: E402
    angle,
    confidence,
    phi_at,
    search_normal,
)

from sfmtool._solve_strips import _SolveStrips  # noqa: E402
from sfmtool._sfmtool import SfmrReconstruction, PatchCloud  # noqa: E402
from sfmtool.sift.file import SiftReader, get_sift_path_from_recon  # noqa: E402


def project(s, i, c):
    pc = s.rot_of[i] @ (c - s.centers[i])
    if pc[2] <= 1e-9:
        return None
    u, v = s.cam_of[i].project(pc[0] / pc[2], pc[1] / pc[2])
    return np.array([u, v])


def sift_offset(s, i, center, u_ax, v_ax, kp, ext):
    """World in-plane offset that moves the patch center so it reprojects onto the
    SIFT keypoint `kp` in view `i` (first-order Jacobian solve). Clamped to ±ext."""
    p0 = project(s, i, center)
    if p0 is None:
        return None, None
    h = ext * 0.05
    ju = (project(s, i, center + h * u_ax) - p0) / h
    jv = (project(s, i, center + h * v_ax) - p0) / h
    J = np.column_stack([ju, jv])
    resid = kp - p0
    try:
        ab = np.linalg.solve(J, resid)
    except np.linalg.LinAlgError:
        return None, None
    wo = ab[0] * u_ax + ab[1] * v_ax
    if np.linalg.norm(wo) > ext:
        wo = wo / np.linalg.norm(wo) * ext
    return wo, float(np.hypot(*resid))


def main() -> None:
    recon_path = sys.argv[1]
    n_points = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    recon = SfmrReconstruction.load(recon_path)
    s = _SolveStrips(recon, recon.workspace_dir, patch=PATCH, extent_factor=EXTENT_FACTOR)
    w = gauss_window(PATCH)
    finite = np.abs(np.asarray(recon.positions_xyzw)[:, 3]) > 1e-9

    # Seed = mean-viewing normal (production's step-1 init before refinement).
    mv = PatchCloud.from_reconstruction(
        recon, normal="mean_viewing", extent_value=EXTENT_FACTOR,
        exclude_points_at_infinity=True,
    )
    seed_of = {int(p): np.asarray(mv[k].normal, float) for k, p in enumerate(mv.point_ids)}

    # Points: finite, well-tracked (need >=4 track views for a meaningful LOO).
    cand = []
    for pid in mv.point_ids:
        pid = int(pid)
        if not finite[pid]:
            continue
        obs = {im: feat for im, feat in s.feat_obs.get(pid, [])}
        if len(obs) >= 4:
            cand.append((pid, obs))
    cand.sort(key=lambda t: len(t[1]), reverse=True)
    pool = cand[:n_points]

    pos_cache: dict[int, np.ndarray] = {}

    def sift_pos(i):
        if i not in pos_cache:
            path = get_sift_path_from_recon(recon, s.names[i])
            pos_cache[i] = np.asarray(SiftReader(str(path)).read_positions(), float)
        return pos_cache[i]

    print(f"{recon_path}")
    print(f"points: {len(pool)} (finite, track>=4)\n")
    hdr = (f"{'pid':>5} {'T':>3} {'resid':>6} | {'move':>6} | "
           f"{'conf proj>sift':>16} | {'LOO proj>sift':>14} | {'Phi proj>sift':>14}")
    print(hdr)
    print("-" * len(hdr))
    moves, dloo, dconf, resids = [], [], [], []
    for pid, obs in pool:
        center = s.positions[pid]
        seed = seed_of[pid]
        T = sorted(obs)
        up = s.rot_of[T[0]].T @ np.array([0.0, -1.0, 0.0])
        ext = s._half(pid)
        from exp_view_localization import patch_frame

        _, u_ax, v_ax, _ = patch_frame(s, pid, seed, T, ext, PATCH)
        c_proj, c_sift, rlist = [], [], []
        for im in T:
            kp = sift_pos(im)[obs[im]]
            wo, r = sift_offset(s, im, center, u_ax, v_ax, kp, ext)
            c_proj.append(center)
            c_sift.append(center if wo is None else center + wo)
            if r is not None:
                rlist.append(r)
        med_resid = float(np.median(rlist)) if rlist else float("nan")

        n_proj, phi_proj = search_normal(s, T, c_proj, seed, up, ext, w)
        n_sift, phi_sift = search_normal(s, T, c_sift, seed, up, ext, w)
        move = angle(n_proj, n_sift)
        conf_proj = confidence(s, T, c_proj, n_proj, up, ext, w)
        conf_sift = confidence(s, T, c_sift, n_sift, up, ext, w)
        _, cores_proj = phi_at(s, T, c_proj, n_proj, up, ext, w)
        _, cores_sift = phi_at(s, T, c_sift, n_sift, up, ext, w)
        loo_proj = loo_zncc(cores_proj, w)
        loo_sift = loo_zncc(cores_sift, w)

        moves.append(move)
        dloo.append(loo_sift - loo_proj)
        dconf.append(conf_sift - conf_proj)
        resids.append(med_resid)
        print(
            f"{pid:>5} {len(T):>3} {med_resid:>6.2f} | {move:>5.1f}° | "
            f"{conf_proj:>7.3f}>{conf_sift:>7.3f} | {loo_proj:>6.3f}>{loo_sift:>6.3f} | "
            f"{phi_proj:>6.3f}>{phi_sift:>6.3f}"
        )

    print("-" * len(hdr))
    print(f"median reprojection residual: {np.median(resids):.2f} img px")
    print(f"median normal move (proj->sift): {np.median(moves):.1f}°")
    print(f"mean ΔLOO (sift-proj): {np.mean(dloo):+.3f}    mean Δconf (sift-proj): {np.mean(dconf):+.3f}")


if __name__ == "__main__":
    main()
