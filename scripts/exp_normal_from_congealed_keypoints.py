#!/usr/bin/env python3
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment: does feeding keypoint positions into normal refinement change the
normal?

Production refines the normal over each view's RAW projection (one shared 3D
center, no per-view offset). `localize_keypoints` later finds a per-view
sub-pixel shift, but that never feeds back into the normal. This probe closes
the loop: it runs the SAME coarse-to-fine normal search twice per point —

  * raw:  every view rendered at the projected center (= production)
  * cong: every view rendered at its keypoint-congealed in-plane offset

— and reports how far the normal moves between them, the Φ-peak sharpness
(confidence), and the leave-one-out ZNCC at each optimum. Same search code both
ways, so the only difference is whether the keypoints are used.

Note: the normal search here is a Python proxy of the Rust kernel — it does not
replicate Rust's frozen-support masking, so its argmax carries a few degrees of
noise vs. the production refiner (see the validation section of
reports/exp/2026-06-21-mvs-normal-refinement.md). The robust signal is the
*paired* ΔLOO, measured at each optimum, not the per-point normal move.

Usage::

    pixi run python scripts/exp_normal_from_congealed_keypoints.py RECON.sfmr [N]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _normal_strip_lib import (  # noqa: E402
    ANG_RANGE,
    EXTENT_FACTOR,
    INIT_STEPS,
    PATCH,
    _phi,
    gauss_window,
    geometric_views,
    irls_weights,
    znorm_stack,
)
from exp_reference_refine import loo_zncc, wpp_to_src  # noqa: E402
from exp_view_localization import (  # noqa: E402
    loo_reference,
    patch_frame,
    render_tile_at_offset,
    search_shift,
)

from sfmtool._solve_strips import _SolveStrips  # noqa: E402
from sfmtool._sfmtool import OrientedPatch, SfmrReconstruction  # noqa: E402
from sfmtool._sfmtool.flow import WarpMap  # noqa: E402

RENDER_RES = 64    # context tile for congealing (needs slide room)
OFF = (RENDER_RES - PATCH) // 2
SEARCH = 6
CONGEAL_ITERS = 5
STEPS = 7          # normal-search grid per axis (matches Rust init_steps)
LEVELS = 3         # coarse-to-fine passes
SHRINK = 0.4


def tangent_basis(n):
    n = n / np.linalg.norm(n)
    a = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = a - n * (a @ n)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    return u, v


def exp_map(n, u, v, a, b):
    """Geodesic step on the unit sphere: tilt n by (a,b) radians in the (u,v) plane."""
    d = a * u + b * v
    th = np.linalg.norm(d)
    if th < 1e-12:
        return n / np.linalg.norm(n)
    return np.cos(th) * n + np.sin(th) * (d / th)


def render_core(s, i, c, n, up, ext):
    """Core PATCH tile for view i, rendered directly at native PATCH resolution:
    shared (n, up) frame, center translated to c (the per-view keypoint offset)."""
    patch = OrientedPatch.from_center_normal(c.tolist(), n.tolist(), up.tolist(), [ext, ext])
    wm = WarpMap.from_patch(patch, s.cam_of[i], s.pose_of[i], PATCH)
    return np.asarray(wm.remap_bilinear(s.image(i)), np.float32)


def phi_at(s, G, centers, n, up, ext, w):
    cores = [render_core(s, G[k], centers[k], n, up, ext) for k in range(len(G))]
    rows = znorm_stack(cores, w)
    if rows.shape[0] < 2:
        return float("nan"), cores
    return _phi(rows, irls_weights(rows)), cores


def search_normal(s, G, centers, seed, up, ext, w):
    """Coarse-to-fine normal search maximizing consensus Phi over the stack."""
    n = seed / np.linalg.norm(seed)
    rng = np.radians(ANG_RANGE)
    best_n, best_phi = n, -np.inf
    for _ in range(LEVELS):
        u, v = tangent_basis(n)
        grid = np.linspace(-rng, rng, STEPS)
        lvl_n, lvl_phi = n, -np.inf
        for a in grid:
            for b in grid:
                nn = exp_map(n, u, v, a, b)
                phi, _ = phi_at(s, G, centers, nn, up, ext, w)
                if phi > lvl_phi:
                    lvl_phi, lvl_n = phi, nn
        n, best_n, best_phi = lvl_n, lvl_n, lvl_phi
        rng *= SHRINK
    return best_n, best_phi


def confidence(s, G, centers, n, up, ext, w):
    """Phi-peak curvature at the optimum (mean of the two tangent second-derivatives);
    higher = more sharply determined normal."""
    u, v = tangent_basis(n)
    h = np.radians(3.0)
    p0, _ = phi_at(s, G, centers, n, up, ext, w)
    pa1, _ = phi_at(s, G, centers, exp_map(n, u, v, h, 0), up, ext, w)
    pa2, _ = phi_at(s, G, centers, exp_map(n, u, v, -h, 0), up, ext, w)
    pb1, _ = phi_at(s, G, centers, exp_map(n, u, v, 0, h), up, ext, w)
    pb2, _ = phi_at(s, G, centers, exp_map(n, u, v, 0, -h), up, ext, w)
    d2 = ((pa1 - 2 * p0 + pa2) + (pb1 - 2 * p0 + pb2)) / (2 * h * h)
    return -d2  # peak => negative curvature => positive confidence


def angle(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    return float(np.degrees(np.arccos(np.clip(abs(a @ b), -1.0, 1.0))))


def congeal_offsets(s, pid, seed, G, ext, w):
    """Per-view world-space in-plane offsets from congealing at the seed normal."""
    center, u_ax, v_ax, wpp = patch_frame(s, pid, seed, G, ext, RENDER_RES)
    nv = len(G)
    acc = np.zeros((nv, 2))
    for _ in range(CONGEAL_ITERS):
        tiles = [
            render_tile_at_offset(s, G[v], center, seed, u_ax, v_ax, acc[v, 0], acc[v, 1], wpp, ext, RENDER_RES)
            for v in range(nv)
        ]
        cores = [t[OFF : OFF + PATCH, OFF : OFF + PATCH] for t in tiles]
        rows = znorm_stack(cores, w)
        deltas = np.zeros((nv, 2))
        for v in range(nv):
            ref = loo_reference(cores, rows, v)
            dx, dy, _, _ = search_shift(tiles[v], ref, OFF, w, SEARCH, True)
            deltas[v] = (dx, dy)
        acc = np.clip(acc + deltas, -SEARCH, SEARCH)
        if float(np.hypot(deltas[:, 0], deltas[:, 1]).mean()) < 0.05:
            break
    wo = [(acc[v, 0] * wpp) * u_ax + (acc[v, 1] * wpp) * v_ax for v in range(nv)]
    src_per_px = wpp_to_src(s, G, center, u_ax, v_ax, wpp)
    kp_src = float(np.median(np.hypot(acc[:, 0], acc[:, 1])) * src_per_px)
    return center, wo, kp_src


def main() -> None:
    recon_path = sys.argv[1]
    n_points = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    recon = SfmrReconstruction.load(recon_path)
    s = _SolveStrips(recon, recon.workspace_dir, patch=PATCH, extent_factor=EXTENT_FACTOR)
    imgs = [s.image(i) for i in range(len(s.names))]
    w = gauss_window(PATCH)
    cidx = {int(p): k for k, p in enumerate(s.cloud.point_ids)}
    finite = np.abs(np.asarray(recon.positions_xyzw)[:, 3]) > 1e-9

    cand = []
    for pid in s.cloud.point_ids:
        pid = int(pid)
        if not finite[pid]:
            continue
        tracked = sorted(set(s.obs.get(pid, [])))
        if len(tracked) < 3:
            continue
        vis = geometric_views(s, pid)
        if len(set(vis) - set(tracked)):
            cand.append((pid, len(set(vis) - set(tracked))))
    cand.sort(key=lambda t: t[1], reverse=True)
    ids = [pid for pid, _ in cand[:n_points]]

    # Track-refine → seed normals (= production's normal). Capture before any
    # further refine mutates the cloud.
    s.cloud.refine_normals(
        recon, imgs, point_ids=ids, resolution=PATCH,
        angular_range_deg=ANG_RANGE, init_steps=INIT_STEPS,
    )
    seed_of = {pid: np.asarray(s.cloud[cidx[pid]].normal, float) for pid in ids}
    sels = {
        int(d["point_id"]): [int(i) for i in np.asarray(d["admitted"]).tolist()]
        for d in s.cloud.select_views(
            recon, imgs, min_relative_zncc=0.7, resolution=PATCH, point_ids=ids
        )
    }
    # Production keypoint localization: kept views per point (occluders dropped).
    # The experiment registers THIS pruned stack — the actual production keypoints.
    kept = {
        int(d["point_id"]): [int(i) for i in np.asarray(d["views"]).tolist()]
        for d in s.cloud.localize_keypoints(
            recon, imgs, view_sets=sels, max_iters=CONGEAL_ITERS, search=float(SEARCH),
            max_shift_px=3.0, min_relative_zncc=0.7, resolution=PATCH,
        )
    }

    # Sanity: a Rust refine over the SAME view set G (zero offsets). My Python
    # raw-search should reproduce this — if it does, the move below is trustworthy.
    sel_set = set(ids)
    track_of = {int(p): sorted(set(s.obs.get(int(p), []))) for p in s.cloud.point_ids}
    vi = []
    for p in s.cloud.point_ids:
        p = int(p)
        if p in sel_set and len(kept.get(p, [])) >= 2:
            vi.append([int(x) for x in sorted(kept[p])])
        else:
            t = track_of[p]
            vi.append([int(x) for x in t] if t else [0])
    rg = s.cloud.refine_normals(
        recon, imgs, point_ids=ids, view_indices=vi, resolution=PATCH,
        angular_range_deg=ANG_RANGE, init_steps=INIT_STEPS,
    )
    rust_G = {pid: np.asarray(rg["normal"])[cidx[pid]] for pid in ids}

    print(f"{recon_path}")
    print(f"points: {len(ids)}\n")
    hdr = (f"{'pid':>5} {'G':>3} {'kp(src)':>7} | {'move':>6} {'sane':>5} | "
           f"{'conf raw>cong':>16} | {'LOO raw>cong':>14} | {'Phi raw>cong':>13}")
    print(hdr)
    print("-" * len(hdr))
    moves, dconf, dloo, kpix, sanes = [], [], [], [], []
    for pid in ids:
        G = sorted(kept[pid])
        if len(G) < 3:
            continue
        seed = seed_of[pid]
        up = s.rot_of[G[0]].T @ np.array([0.0, -1.0, 0.0])
        ext_core = s._half(pid)
        ext_ctx = s._half(pid) * (RENDER_RES / PATCH)
        center, wo, kp_src = congeal_offsets(s, pid, seed, G, ext_ctx, w)
        c_raw = [center for _ in G]
        c_cong = [center + o for o in wo]

        n_raw, phi_raw = search_normal(s, G, c_raw, seed, up, ext_core, w)
        n_cong, phi_cong = search_normal(s, G, c_cong, seed, up, ext_core, w)

        move = angle(n_raw, n_cong)
        sane = angle(n_raw, rust_G[pid])  # Python raw-search vs Rust-over-G
        conf_raw = confidence(s, G, c_raw, n_raw, up, ext_core, w)
        conf_cong = confidence(s, G, c_cong, n_cong, up, ext_core, w)
        _, raw_cores = phi_at(s, G, c_raw, n_raw, up, ext_core, w)
        _, cong_cores = phi_at(s, G, c_cong, n_cong, up, ext_core, w)
        loo_raw = loo_zncc(raw_cores, w)
        loo_cong = loo_zncc(cong_cores, w)

        moves.append(move)
        sanes.append(sane)
        dconf.append(conf_cong - conf_raw)
        dloo.append(loo_cong - loo_raw)
        kpix.append(kp_src)
        print(
            f"{pid:>5} {len(G):>3} {kp_src:>7.2f} | {move:>5.1f}° {sane:>4.1f}° | "
            f"{conf_raw:>7.3f}>{conf_cong:>7.3f} | {loo_raw:>6.3f}>{loo_cong:>6.3f} | "
            f"{phi_raw:>6.3f}>{phi_cong:>6.3f}"
        )

    print("-" * len(hdr))
    print(f"median keypoint shift: {np.median(kpix):.2f} src px   "
          f"(search sanity vs Rust-over-G: median {np.median(sanes):.1f}°)")
    print(f"median normal move (raw->cong): {np.median(moves):.1f}°")
    print(f"mean ΔLOO (cong-raw): {np.mean(dloo):+.3f}    mean Δconf (cong-raw): {np.mean(dconf):+.3f}")


if __name__ == "__main__":
    main()
