# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment: rapid pinhole-camera estimate from cluster patches.

The first stage of a divide-and-conquer bootstrap: from a workspace holding
a `*-clusters-patches.matches` file, get to a good shared SIMPLE_PINHOLE
estimate (focal + a small set of posed views) as fast as possible.  This is
`exp_pinhole_bootstrap.py` cut down to the part before full growth:

  1. Load patch clusters (span >= 2), keep only the best SFMTOOL_MAX_CL
     clusters by span (the wide clusters carry both the covisibility signal
     and the focal observability).
  2. Covisibility seed groups -> affine ALS factorization + Tomasi-Kanade
     metric upgrade of the best two candidate windows.
  3. Fixed-focal scan over a small grid: seed a perspective solve per
     candidate (both reflection hypotheses), grow by P3P resection to a
     small image cap (SFMTOOL_SCAN_CAP), rank by inlier fraction.
  4. Release f in a staged BA on the winner's subset -> final estimate.

Run: pixi run -e dev python scripts/exp_fast_pinhole.py <workspace> [ref.sfmr]

Prints the focal + camera errors vs the reference solve (when one exists)
and writes `<workspace>/fast-pinhole.json` with the estimate for later
stages to consume.
"""

import itertools
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import scipy.optimize
import scipy.sparse
from scipy.spatial.transform import Rotation

from sfmtool._sfmtool.geometry import (
    CameraIntrinsics,
    estimate_absolute_pose,
    factorize_affine,
    inlier_fraction as _inlier_fraction,
    refine_absolute_pose as _refine_absolute_pose,
    reprojection_residuals as _reprojection_residuals,
)

WS = Path(sys.argv[1] if len(sys.argv) > 1 else "e_seoul_ws")
REF = Path(sys.argv[2]) if len(sys.argv) > 2 else None
_T0 = time.perf_counter()

# BA budget: only the best MAX_CL clusters (highest span first) enter the
# bundle adjustments.  Covisibility, factorization, resection, and
# triangulation always see every usable cluster — a capped working set
# starves the seed window (dino's high-span prefix is junk-dominated and a
# 2000-cluster load cap collapsed its factorization outright).
MAX_CL = int(os.environ.get("SFMTOOL_MAX_CL", "3000"))
SCAN_CAP = int(os.environ.get("SFMTOOL_SCAN_CAP", "8"))  # core images per try
# Wide-baseline shell: after the covisibility-driven core, resect up to this
# many far images (widest viewpoint angles).  The focal is unobservable on a
# near-affine core — a sliver of a long orbit fits ANY focal at high inlier
# fraction (living-room: 86% inliers at f 45% off) — and the shell is what
# makes the scan discriminate and the released f converge.
SHELL = int(os.environ.get("SFMTOOL_SHELL", "6"))
# Per-image BA row cap: each image keeps only its best observations (by
# admission rank), so BA cost stays flat as the shell widens the image set.
OBS_PER_IMG = int(os.environ.get("SFMTOOL_OBS_PER_IMG", "250"))
# Anchor-cluster budget for the photometric verification pass.
N_ANCHORS = int(os.environ.get("SFMTOOL_ANCHORS", "400"))
F_GRID = [0.55, 0.7, 0.9, 1.2, 1.6]  # focal candidates, units of max(w, h)

# Canonical camera frame (-Z forward, +Y up) throughout; full-pixel
# observations; shared SIMPLE_PINHOLE with the principal point at the centre.
_CAM_WH = None


def elapsed():
    return time.perf_counter() - _T0


def make_cam(f):
    w, h = _CAM_WH
    return CameraIntrinsics.from_dict(
        {
            "model": "SIMPLE_PINHOLE",
            "width": int(w),
            "height": int(h),
            "parameters": {
                "focal_length": float(f),
                "principal_point_x": w / 2.0,
                "principal_point_y": h / 2.0,
            },
        }
    )


# ── Data loading ─────────────────────────────────────────────────────────────


def load_clusters():
    """Flat observation arrays of the best MAX_CL clusters by span.

    Refined member positions come from the stored affines' last column
    (absolute keypoint position, .matches format v4+); no per-image reads.
    """
    from sfmtool._sfmtool.io import read_matches

    override = os.environ.get("SFMTOOL_MATCHES")
    patches = (
        [Path(override)]
        if override
        else sorted(WS.glob("matches/*-clusters-patches.matches"))
    )
    data = read_matches(patches[0])
    names = list(data["image_names"])
    dims = [(int(w), int(h)) for w, h in np.asarray(data["image_dims"])]

    starts = np.asarray(data["cluster_starts"])
    mi = np.asarray(data["member_images"])
    mf = np.asarray(data["member_features"])
    st = np.asarray(data["member_status"])
    refs = np.asarray(data["reference_members"])
    aff = np.asarray(data["member_affines"])
    cons = np.asarray(data["member_consistency_residual"], dtype=np.float64)

    # (span, cluster id, member selection, worst finite warp-consistency
    # residual — lower = better; inf where no member entered the fit)
    usable = []
    for c in range(len(starts) - 1):
        lo, hi = int(starts[c]), int(starts[c + 1])
        if refs[c] == np.iinfo(np.uint32).max:
            continue
        sel = np.nonzero((st[lo:hi] == 0) | (st[lo:hi] == 1))[0] + lo
        span = len(np.unique(mi[sel]))
        if span >= 2:
            cq = cons[sel]
            cq = cq[np.isfinite(cq)]
            usable.append((span, c, sel, float(cq.max()) if len(cq) else np.inf))
    spans = np.array([t[0] for t in usable])
    cids = np.array([t[1] for t in usable])
    order = np.lexsort((cids, -spans))
    rank = np.empty(len(usable), np.int64)
    rank[order] = np.arange(len(usable))

    obs_c, obs_i, obs_f, obs_uv, adm_rank, quality = [], [], [], [], [], []
    for n_cl, k in enumerate(sorted(range(len(usable)), key=lambda k: cids[k])):
        _span, _c, sel, q = usable[k]
        adm_rank.append(rank[k])
        quality.append(q)
        for m in sel:
            obs_c.append(n_cl)
            obs_i.append(int(mi[m]))
            obs_f.append(int(mf[m]))
            obs_uv.append(aff[m, :, 2])
    return {
        "names": names,
        "dims": dims,
        "obs_c": np.asarray(obs_c),
        "obs_i": np.asarray(obs_i),
        "obs_f": np.asarray(obs_f),
        "obs_uv": np.asarray(obs_uv, dtype=np.float64),
        "adm_rank": np.asarray(adm_rank, dtype=np.int64),
        "cl_quality": np.asarray(quality, dtype=np.float64),
        "n_img": len(names),
        "n_cl": len(usable),
    }


# ── Seed: covisibility grouping + affine factorization ───────────────────────


def build_covisibility(obs_c, obs_i, n_img, n_cl):
    from sfmtool._sfmtool.matching import ClusterCovisibility

    starts = np.searchsorted(obs_c, np.arange(n_cl + 1)).astype(np.uint32)
    return ClusterCovisibility.from_arrays(starts, obs_i.astype(np.uint32), n_img)


def window_spans(obs_c, obs_i, imgs, min_span):
    """Observation selection for clusters seen in >= min_span window images."""
    inw = np.isin(obs_i, imgs)
    cl, il = obs_c[inw], np.searchsorted(imgs, obs_i[inw])
    span = np.zeros(cl.max() + 1 if len(cl) else 1, int)
    for c in np.unique(cl):
        span[c] = len(np.unique(il[cl == c]))
    sel = inw.copy()
    sel[inw] = np.isin(cl, np.nonzero(span >= min_span)[0])
    uniq, c2 = np.unique(obs_c[sel], return_inverse=True)
    return sel, np.searchsorted(imgs, obs_i[sel]), uniq, c2


def factorize_window(obs_c, obs_i, u, imgs, min_span=3):
    """ALS affine factorization + metric upgrade of a candidate window.

    Returns (metric hypotheses as (rot, scale, t_aff), used mask, span-2
    selection for the window mini-BA) or None when too sparse.
    """
    sel, il, uniq, c2 = window_spans(obs_c, obs_i, imgs, min_span)
    if sel.sum() < 30:
        return None
    fac = factorize_affine(
        c2.astype(np.uint32),
        il.astype(np.uint32),
        np.ascontiguousarray(u[sel]),
        len(imgs),
        len(uniq),
    )
    used = np.asarray(fac.used_images)
    t_aff = np.asarray(fac.translations)
    upgraded = fac.metric_upgrade()
    hyps = []
    if upgraded is not None:
        for hyp in upgraded:
            rot = np.asarray(hyp.rotations)
            scale = np.asarray(hyp.scales)
            if (scale[used] > 0).all():
                hyps.append((rot, scale, t_aff))
    return hyps, used, window_spans(obs_c, obs_i, imgs, 2)


def perspective_init(rot, scale, t_cam, used, f0):
    """Weak-perspective -> pinhole poses in the canonical camera frame."""
    rot_can = rot * np.array([1.0, -1.0, -1.0])[None, :, None]
    trans = np.zeros((len(rot), 3))
    for i in np.nonzero(used)[0]:
        trans[i] = [t_cam[i, 0] / scale[i], -t_cam[i, 1] / scale[i], -f0 / scale[i]]
    return rot_can, trans


# ── Geometry kernels ─────────────────────────────────────────────────────────


def triangulate(obs_c, obs_i, u, rot, trans, used, n_cl, f):
    """Ray-midpoint triangulation from the posed images (< 2 views: NaN)."""
    from sfmtool._sfmtool.analysis import triangulate_batch

    pts = np.full((n_cl, 3), np.nan)
    sel = used[obs_i]
    if not sel.any():
        return pts
    oc, oi, uv = obs_c[sel], obs_i[sel], u[sel]
    d_loc = make_cam(f).pixel_to_ray_batch(np.ascontiguousarray(uv))
    dirs = np.einsum("nji,nj->ni", rot[oi], d_loc)
    centers = -np.einsum("nji,nj->ni", rot[oi], trans[oi])
    uniq, counts = np.unique(oc, return_counts=True)
    offsets = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
    result = triangulate_batch(
        np.ascontiguousarray(dirs), np.ascontiguousarray(centers), offsets
    )
    good = counts >= 2
    pts[uniq[good]] = np.asarray(result["points"])[good]
    return pts


def p3p_resect(uv, x_pts, f0):
    """RANSAC P3P absolute pose; returns (rvec, tvec, inlier mask) or None."""
    ans = estimate_absolute_pose(
        np.ascontiguousarray(uv),
        np.ascontiguousarray(x_pts),
        camera=make_cam(f0),
        max_error_px=4.0,
        seed=0,
    )
    if ans is None:
        return None
    q = np.asarray(ans["quaternion_wxyz"])
    rv = Rotation.from_quat(q[[1, 2, 3, 0]]).as_rotvec()
    return rv, np.asarray(ans["translation"]), np.asarray(ans["inliers"], dtype=bool)


def pose_refine(uv, x_pts, rv0, tv0, f):
    """Trimmed pose-only refinement (native)."""
    q0 = Rotation.from_rotvec(rv0).as_quat()[[3, 0, 1, 2]]
    out = _refine_absolute_pose(
        make_cam(f),
        np.ascontiguousarray(uv, dtype=np.float64),
        np.ascontiguousarray(x_pts, dtype=np.float64),
        q0,
        np.ascontiguousarray(tv0, dtype=np.float64),
        5,
        0.6,
        3.0,
    )
    q = np.asarray(out["quaternion_wxyz"])
    rv = Rotation.from_quat(q[[1, 2, 3, 0]]).as_rotvec()
    return rv, np.asarray(out["translation"]), float(out["inlier_fraction"])


def reproj_res_one(cam, rvec_i, tvec_i, x_pts, uv):
    q = Rotation.from_rotvec(rvec_i).as_quat()[[3, 0, 1, 2]][None, :]
    n = len(uv)
    return _reprojection_residuals(
        cam,
        q,
        np.ascontiguousarray(tvec_i, dtype=np.float64)[None, :],
        np.ascontiguousarray(x_pts, dtype=np.float64),
        np.ascontiguousarray(uv, dtype=np.float64),
        np.zeros(n, np.uint32),
        np.arange(n, dtype=np.uint32),
        1e6,
    )


# ── Bundle adjustment ────────────────────────────────────────────────────────


def project(cam, rvec, tvec, pts, obs_i, obs_c):
    xc = Rotation.from_rotvec(rvec[obs_i]).apply(pts[obs_c]) + tvec[obs_i]
    return cam.ray_to_pixel_batch(np.ascontiguousarray(xc)), -xc[:, 2]


def solve_round(obs_c, obs_i, u, f, rvec, tvec, pts, f_scale, opt_f, max_nfev):
    """One robust sparse least-squares pass (compact-indexed)."""
    uniq, oc2 = np.unique(obs_c, return_inverse=True)
    uimg, oi2 = np.unique(obs_i, return_inverse=True)
    n_pt, n_im = len(uniq), len(uimg)
    nf = 1 if opt_f else 0

    def unpack(x):
        fv = x[0] if opt_f else f
        rv = x[nf : nf + 3 * n_im].reshape(-1, 3)
        tv = x[nf + 3 * n_im : nf + 6 * n_im].reshape(-1, 3)
        pv = x[nf + 6 * n_im :].reshape(-1, 3)
        return fv, rv, tv, pv

    base_cam = None if opt_f else make_cam(f)
    oi2_u32, oc2_u32 = oi2.astype(np.uint32), oc2.astype(np.uint32)
    u_c = np.ascontiguousarray(u, dtype=np.float64)

    def resid(x):
        fv, rv, tv, pv = unpack(x)
        cam = base_cam if base_cam is not None else make_cam(fv)
        q = Rotation.from_rotvec(rv).as_quat()[:, [3, 0, 1, 2]]
        return _reprojection_residuals(
            cam,
            np.ascontiguousarray(q),
            np.ascontiguousarray(tv),
            np.ascontiguousarray(pv),
            u_c,
            oi2_u32,
            oc2_u32,
            1e6,
        ).ravel()

    n_obs = len(oc2)
    spar = scipy.sparse.lil_matrix(
        (2 * n_obs, nf + 6 * n_im + 3 * n_pt), dtype=np.uint8
    )
    rows = np.arange(2 * n_obs)
    if opt_f:
        spar[:, 0] = 1
    for k in range(3):
        spar[rows, nf + 3 * np.repeat(oi2, 2) + k] = 1
        spar[rows, nf + 3 * n_im + 3 * np.repeat(oi2, 2) + k] = 1
        spar[rows, nf + 6 * n_im + 3 * np.repeat(oc2, 2) + k] = 1

    x0 = np.concatenate(
        [
            [f] if opt_f else [],
            rvec[uimg].ravel(),
            tvec[uimg].ravel(),
            pts[uniq].ravel(),
        ]
    )
    sol = scipy.optimize.least_squares(
        resid,
        x0,
        jac_sparsity=spar,
        loss="soft_l1",
        f_scale=f_scale,
        max_nfev=max_nfev,
        x_scale="jac",
        verbose=0,
    )
    f, rv, tv, p_new = unpack(sol.x)
    rvec, tvec, out = rvec.copy(), tvec.copy(), pts.copy()
    rvec[uimg], tvec[uimg], out[uniq] = rv, tv, p_new
    return f, rvec, tvec, out


def bundle_adjust(
    obs_c,
    obs_i,
    u,
    rvec,
    tvec,
    pts,
    f0,
    n_img,
    n_cl,
    opt_f,
    schedule=((50.0, 5.0), (12.0, 2.0), (4.0, 1.0)),
    max_nfev=60,
):
    """Staged robust BA: trim, solve, re-triangulate between rounds."""
    f = f0
    all_img = np.ones(n_img, bool)
    for rnd, (thresh, f_scale) in enumerate(schedule):
        if rnd > 0:
            rot_now = Rotation.from_rotvec(rvec).as_matrix()
            pts = triangulate(obs_c, obs_i, u, rot_now, tvec, all_img, n_cl, f)
        proj, depth = project(make_cam(f), rvec, tvec, pts, obs_i, obs_c)
        rn = np.linalg.norm(proj - u, axis=1)
        with np.errstate(invalid="ignore"):
            keep = (rn < thresh) & (depth > 1e-3 * f) & ~np.isnan(pts[obs_c, 0])
        surv = np.bincount(obs_c[keep], minlength=n_cl)
        keep &= surv[obs_c] >= 2
        if keep.sum() < 12:  # degenerate (e.g. a wildly wrong focal)
            return f, rvec, tvec, pts, np.full(len(obs_c), np.inf), 0.0
        f, rvec, tvec, pts = solve_round(
            obs_c[keep],
            obs_i[keep],
            u[keep],
            f,
            rvec,
            tvec,
            pts,
            f_scale,
            opt_f,
            max_nfev,
        )
    proj, depth = project(make_cam(f), rvec, tvec, pts, obs_i, obs_c)
    rn = np.linalg.norm(proj - u, axis=1)
    res = np.where(np.isnan(rn), np.inf, rn)
    return f, rvec, tvec, pts, res, float((res < 2.0).mean())


# ── Growth (scan-sized) ──────────────────────────────────────────────────────

_RANK_O = None  # per-observation admission rank (set in main)


def ba_rows(live, obs_i):
    """Per-image cap on BA rows: keep each image's best OBS_PER_IMG
    observations by admission rank, so BA cost stays flat in image count."""
    idx = np.nonzero(live)[0]
    keep = live.copy()
    for i in np.unique(obs_i[idx]):
        rows = idx[obs_i[idx] == i]
        if len(rows) > OBS_PER_IMG:
            keep[rows[np.argsort(_RANK_O[rows], kind="stable")[OBS_PER_IMG:]]] = (
                False
            )
    return keep


def grow_to_cap(seed, f0, obs_c, obs_i, u, n_img, n_cl, cap, bam):
    """Seed perspective solves and P3P-grow each to ``cap`` images.

    Both reflection hypotheses of the metric upgrade fit a near-affine seed
    window almost equally well, so the top two seed candidates (by the seed
    mini-BA inlier fraction) BOTH grow — the mirror solution falls behind
    once wider-baseline views join.  Yields the grown states for the caller
    to rank.

    Minimal next-best-view loop: no force-accept or retry machinery — an
    image that fails resection or the inlier gate is blocked for good (the
    scan ranks candidates; it does not need completion).
    """
    cands = []
    for imgs, wd in seed:
        if wd is None:
            continue
        hyps, used, (sel, il, cl_ids, c2) = wd
        for rot0, scale, t_aff in hyps:
            rot_can, trans0 = perspective_init(rot0, scale, t_aff, used, f0)
            pts_w = triangulate(c2, il, u[sel], rot_can, trans0, used, len(cl_ids), f0)
            ok = ~np.isnan(pts_w[:, 0])[c2] & used[il]
            if ok.sum() < 30:
                continue
            rvec_w = Rotation.from_matrix(rot_can).as_rotvec()
            _, rvw, tvw, p_w, _, inl = bundle_adjust(
                c2[ok],
                il[ok],
                u[sel][ok],
                rvec_w,
                trans0,
                pts_w,
                f0,
                len(imgs),
                len(cl_ids),
                opt_f=False,
                schedule=((30.0, 3.0), (8.0, 1.5)),
                max_nfev=30,
            )
            cands.append((inl, imgs, used, cl_ids, rvw, tvw, p_w))
    cands.sort(key=lambda t: -t[0])
    return [_grow_one(c, f0, obs_c, obs_i, u, n_img, n_cl, cap, bam) for c in cands[:2]]


def core_parallax(rvec, tvec, pts, posed, obs_c, obs_i):
    """Median over triangulated points of the widest ray angle between the
    posed views observing them, in degrees.

    A covisibility-picked seed can be a zero-baseline segment (a video's
    most-mutually-covisible frames are where the camera moved LEAST —
    DinoLedge's seed was a near-static clip at the end of the walk).  Such
    a core fits any focal at high inlier fraction while its depths are
    unusable, so growth and the focal scan need a parallax gate, not a
    reprojection one."""
    valid = ~np.isnan(pts[:, 0])[obs_c] & posed[obs_i]
    if not valid.any():
        return 0.0
    oc, oi = obs_c[valid], obs_i[valid]
    rot = Rotation.from_rotvec(rvec).as_matrix()
    centers = -np.einsum("nji,nj->ni", rot[oi], tvec[oi])
    d = pts[oc] - centers
    d /= np.linalg.norm(d, axis=1, keepdims=True) + 1e-12
    order = np.argsort(oc, kind="stable")
    oc_s, d_s = oc[order], d[order]
    uniq, starts = np.unique(oc_s, return_index=True)
    first = np.repeat(d_s[starts], np.diff(np.append(starts, len(oc_s))), axis=0)
    cosang = np.clip((d_s * first).sum(1), -1.0, 1.0)
    widest = np.minimum.reduceat(cosang, starts)
    return float(np.degrees(np.median(np.arccos(widest))))


def _grow_one(cand, f0, obs_c, obs_i, u, n_img, n_cl, cap, bam):
    _, imgs, used, cl_ids, rvw, tvw, p_w = cand

    rvec = np.zeros((n_img, 3))
    tvec = np.tile([0.0, 0.0, -f0], (n_img, 1))
    posed = np.zeros(n_img, bool)
    pts = np.full((n_cl, 3), np.nan)
    for k, i in enumerate(imgs):
        if used[k]:
            rvec[i], tvec[i], posed[i] = rvw[k], tvw[k], True
    pts[cl_ids] = p_w

    cam0 = make_cam(f0)

    def refill(pts, rvec, tvec):
        # Triangulate clusters that lack a point but have >= 2 posed views
        # (also restores the points a budget-restricted BA wiped to NaN).
        need = np.isnan(pts[:, 0])[obs_c] & posed[obs_i]
        if need.any():
            uniq, c2n = np.unique(obs_c[need], return_inverse=True)
            rot = Rotation.from_rotvec(rvec).as_matrix()
            pts[uniq] = triangulate(
                c2n, obs_i[need], u[need], rot, tvec, posed, len(uniq), f0
            )
        return pts

    accepted, blocked, since_ba = [], set(), 0
    while posed.sum() < cap:
        cand = ~posed[obs_i] & ~np.isnan(pts[obs_c, 0])
        cnt = np.bincount(obs_i[cand], minlength=n_img)
        for j in blocked:
            cnt[j] = 0
        i = int(np.argmax(cnt))
        if cnt[i] < 6:
            break
        s = (obs_i == i) & ~np.isnan(pts[obs_c, 0])
        found = None
        p3p = p3p_resect(u[s], pts[obs_c[s]], f0)
        if p3p is not None and int(p3p[2].sum()) >= 12:
            rv0, tv0, mask = p3p
            rv, tv, _ = pose_refine(u[s][mask], pts[obs_c[s]][mask], rv0, tv0, f0)
            res = reproj_res_one(cam0, rv, tv, pts[obs_c[s]], u[s])
            found = (_inlier_fraction(res, 3.0), rv, tv)
        if found is None or (accepted and found[0] < 0.35 * float(np.median(accepted))):
            blocked.add(i)
            continue
        accepted.append(found[0])
        rvec[i], tvec[i], posed[i] = found[1], found[2], True
        pts = refill(pts, rvec, tvec)
        since_ba += 1
        if since_ba >= 3 and posed.sum() < cap:
            since_ba = 0
            live = ba_rows(posed[obs_i] & ~np.isnan(pts[obs_c, 0]) & bam, obs_i)
            _, rvec, tvec, pts, _, _ = bundle_adjust(
                obs_c[live],
                obs_i[live],
                u[live],
                rvec,
                tvec,
                pts,
                f0,
                n_img,
                n_cl,
                opt_f=False,
                schedule=((30.0, 3.0), (8.0, 1.5)),
                max_nfev=30,
            )
            pts = refill(pts, rvec, tvec)

    med_inl = float(np.median(accepted)) if accepted else 1.0
    return rvec, tvec, pts, posed, med_inl


def widen(rvec, tvec, pts, posed, f0, obs_c, obs_i, u, n_img, n_cl, bam, gate):
    """Ladder-widen a converged fixed-f state for focal observability.

    Structure triangulated from a near-affine core has depth errors too
    large for far views to resect against directly (a 12-image jump probe
    on the console orbit failed at every focal).  Instead each rung resects
    the FARTHEST currently-viable image (weakest covisibility link that
    passes the inlier gate), then re-triangulates and bundle-adjusts — the
    reachable arc grows with every rung, so a handful of rungs spans an
    orbit that incremental most-covisible growth needs dozens of images to
    cross.  The focal stays fixed; the caller releases it afterwards.
    """
    cam0 = make_cam(f0)

    def refill(pts):
        need = np.isnan(pts[:, 0])[obs_c] & posed[obs_i]
        if need.any():
            uniq, c2n = np.unique(obs_c[need], return_inverse=True)
            rot = Rotation.from_rotvec(rvec).as_matrix()
            pts[uniq] = triangulate(
                c2n, obs_i[need], u[need], rot, tvec, posed, len(uniq), f0
            )
        return pts

    rejected = set()
    for _rung in range(SHELL):
        valid = ~np.isnan(pts[obs_c, 0])
        cnt = np.bincount(obs_i[~posed[obs_i] & valid], minlength=n_img)
        for j in rejected:
            cnt[j] = 0
        pool = np.nonzero((cnt >= 30) & ~posed)[0]
        # Weakest link = farthest first, falling back toward nearer views
        # until one resects (the farthest reachable image extends the arc,
        # and the reachable radius grows as rungs accumulate).  Log-spaced
        # sampling keeps the near end of the pool in reach — a rung whose
        # far candidates are all junk-connected must still make progress.
        pool = pool[np.argsort(cnt[pool])]
        if len(pool) > 20:
            pool = pool[np.unique(np.geomspace(1, len(pool), 20).astype(int) - 1)]
        hit = None
        for j in pool:
            s = (obs_i == j) & valid
            p3p = p3p_resect(u[s], pts[obs_c[s]], f0)
            if p3p is None or int(p3p[2].sum()) < 12:
                continue
            rv0, tv0, mask = p3p
            rv, tv, _ = pose_refine(u[s][mask], pts[obs_c[s]][mask], rv0, tv0, f0)
            res = reproj_res_one(cam0, rv, tv, pts[obs_c[s]], u[s])
            # Relative gate (like core growth): an absolute floor rejects
            # every candidate when the probe focal is far from true and the
            # whole fit sits at a low but consistent inlier level.
            if _inlier_fraction(res, 3.0) < gate:
                continue
            hit = (j, rv, tv)
            break
        if hit is None:
            break
        # Verified acceptance: a far P3P against depth-noisy points can
        # find a junk consensus that wrecks the BA (wide span, broken
        # geometry).  Accept, BA, then require the image to have SURVIVED
        # the adjustment; revert and blacklist it otherwise.  (A global
        # "did the old images keep their fit" check does NOT work here: a
        # legitimate far rung worsens the old fit by design — breaking the
        # bas-relief compensation is what it is for.)
        saved = (rvec.copy(), tvec.copy(), pts.copy(), posed.copy())
        j, rvec[j], tvec[j] = hit
        posed[j] = True
        pts = refill(pts)
        live = ba_rows(posed[obs_i] & ~np.isnan(pts[obs_c, 0]) & bam, obs_i)
        _, rvec, tvec, pts, _, _ = bundle_adjust(
            obs_c[live],
            obs_i[live],
            u[live],
            rvec,
            tvec,
            pts,
            f0,
            n_img,
            n_cl,
            opt_f=False,
            schedule=((12.0, 2.0), (4.0, 1.0)),
            max_nfev=30,
        )
        pts = refill(pts)
        s = (obs_i == j) & ~np.isnan(pts[obs_c, 0])
        res = reproj_res_one(cam0, rvec[j], tvec[j], pts[obs_c[s]], u[s])
        if _inlier_fraction(res, 3.0) < gate:
            rvec, tvec, pts, posed = saved
            rejected.add(int(j))
    return rvec, tvec, pts, posed


# ── Photometric verification (embed-patches machinery) ──────────────────────


def localize_anchors(names, sub, rvec, tvec, f0, pts_a, tr_a, tr_img, tr_feat):
    """Congeal-localize anchor keypoints across an image subset.

    Builds an in-memory reconstruction over ``sub`` (image indexes, all
    posed), a feature-scaled mean-viewing-normal patch cloud over the
    anchor points, and an in-memory pyramid set, then localizes every
    anchor in EVERY subset view (a patch tile is registered against the
    leave-one-out consensus of the other views' tiles).  Returns
    (anchor idx, full image idx, keypoint xy) arrays of the KEPT views —
    appearance-verified observations.  Cost: ~1-2 s for 400 anchors over
    ~15 4K frames, pyramids included.
    """
    import cv2

    from sfmtool._sfmtool import ImagePyramidSet, PatchCloud, SfmrReconstruction
    from sfmtool._workspace import load_workspace_config
    from sfmtool.colmap.io import (
        _build_sfmr_data_dict,
        _resolve_workspace_and_sift,
        build_metadata,
        finite_positions_xyzw,
    )

    sub_names = [names[int(g)] for g in sub]
    q = Rotation.from_rotvec(rvec[sub]).as_quat()[:, [3, 0, 1, 2]]
    counts = np.bincount(tr_a, minlength=len(pts_a)).astype(np.uint32)
    wsdir, _contents, resolved, ft_hashes, sc_hashes, thumbs = (
        _resolve_workspace_and_sift(sub_names, WS.resolve())
    )
    metadata = build_metadata(
        workspace_dir=wsdir,
        output_path=WS.resolve() / "sfmr" / "fast-pinhole-anchors.sfmr",
        workspace_config=load_workspace_config(wsdir),
        operation="fast_pinhole_anchors",
        tool_name="sfmtool",
        tool_options={},
        image_count=len(sub_names),
        point_count=len(pts_a),
        observation_count=int(counts.sum()),
        camera_count=1,
    )
    sfmr_dict = _build_sfmr_data_dict(
        cameras=[make_cam(f0)],
        image_names=resolved,
        camera_indexes=np.zeros(len(sub_names), dtype=np.uint32),
        quaternions_wxyz=q,
        translations_xyz=tvec[sub],
        positions_xyzw=finite_positions_xyzw(pts_a),
        colors_rgb=np.zeros((len(pts_a), 3), np.uint8),
        reprojection_errors=np.zeros(len(pts_a), np.float32),
        track_image_indexes=tr_img.astype(np.uint32),
        track_feature_indexes=tr_feat.astype(np.uint32),
        point_indexes=tr_a.astype(np.uint32),
        observation_counts=counts,
        feature_tool_hashes=ft_hashes,
        sift_content_hashes=sc_hashes,
        thumbnails=thumbs,
        metadata=metadata,
    )
    recon = SfmrReconstruction.from_data(wsdir, sfmr_dict)
    imgs = [
        np.ascontiguousarray(cv2.imread(str(wsdir / n), cv2.IMREAD_COLOR))
        for n in sub_names
    ]
    pyrset = ImagePyramidSet(recon, imgs)
    # Patch size must track the FEATURE scale: a fixed pixel radius that is
    # fine on a 480 px frame is hopelessly small on 4K (a 12 px tile has no
    # discriminative texture and a patch-grid search budget of a few source
    # px, so localization can neither reach nor reject anything).
    cloud = PatchCloud.from_reconstruction(
        recon, normal="mean_viewing", extent="feature_size", extent_value=2.5
    )
    all_views = list(range(len(sub_names)))
    results = cloud.localize_keypoints(
        recon,
        pyrset,
        view_sets=dict.fromkeys(range(len(pts_a)), all_views),
        max_shift_px=60.0,
        search=12.0,
        min_relative_zncc=0.6,
    )
    a_idx, i_idx, uv = [], [], []
    for r in results:
        pid = int(r["point_index"])
        views = np.asarray(r["views"])
        kps = np.asarray(r["keypoints"], dtype=np.float64)
        for k, v in enumerate(views):
            a_idx.append(pid)
            i_idx.append(int(sub[int(v)]))
            uv.append(kps[k])
    return (
        np.asarray(a_idx),
        np.asarray(i_idx),
        np.asarray(uv, dtype=np.float64).reshape(-1, 2),
    )


# ── Evaluation ───────────────────────────────────────────────────────────────


def compare_to_reference(names, rvec, tvec, f_est, mask):
    """Camera errors vs the first non-bootstrap solve in the workspace.

    A posed SUBSET can have nearly-degenerate camera centers (a short arc of
    a long orbit), which leaves the center-based similarity alignment a free
    rotation about the arc — so the gauge rotation for the ROTATION errors
    is fitted from the camera rotations (well-conditioned always), and the
    center errors use the free similarity (its own best case).
    """
    names = [n for j, n in enumerate(names) if mask[j]]
    rvec, tvec = rvec[mask], tvec[mask]
    ref_files = (
        [REF]
        if REF is not None
        else sorted(
            p
            for p in WS.glob("sfmr/*.sfmr")
            if "bootstrap" not in p.name and "fast-pinhole" not in p.name
        )
    )
    if not ref_files:
        print("no reference solve found; skipping comparison")
        return
    from sfmtool._sfmtool import SfmrReconstruction
    from sfmtool._sfmtool.analysis import estimate_alignment_rs

    ref = SfmrReconstruction.load(ref_files[0])
    ref_names = list(ref.image_names)
    common = [n for n in names if n in ref_names]
    if len(common) < 3:
        print(f"only {len(common)} common images with {ref_files[0].name}; skipping")
        return

    def centers_rots(qs, ts, order):
        rs = Rotation.from_quat(np.asarray(qs)[order][:, [1, 2, 3, 0]]).as_matrix()
        return -np.einsum("nij,ni->nj", rs, np.asarray(ts)[order]), rs

    q_wxyz = Rotation.from_rotvec(rvec).as_quat()[:, [3, 0, 1, 2]]
    ei = np.array([names.index(n) for n in common])
    ri = np.array([ref_names.index(n) for n in common])
    c_est, r_est = centers_rots(q_wxyz, tvec, ei)
    c_ref, r_ref = centers_rots(ref.quaternions_wxyz, ref.translations, ri)

    # Gauge rotation from the rotations: argmin_g sum ||R_est_i g - R_ref_i||.
    u_svd, _s, vt = np.linalg.svd(np.einsum("nji,njk->ik", r_est, r_ref))
    if np.linalg.det(u_svd @ vt) < 0:
        u_svd[:, 2] *= -1.0
    g = u_svd @ vt
    rot_err = Rotation.from_matrix(
        np.einsum("nij,nkj->nik", r_ref, np.einsum("nij,jk->nik", r_est, g))
    ).magnitude() * (180 / np.pi)

    tf = estimate_alignment_rs(
        np.ascontiguousarray(c_est, dtype=np.float64),
        np.ascontiguousarray(c_ref, dtype=np.float64),
    )
    c_fit = tf.apply_to_points(np.ascontiguousarray(c_est, dtype=np.float64))
    diam = np.max(np.linalg.norm(c_ref[:, None, :] - c_ref[None, :, :], axis=2))
    cen_err = np.linalg.norm(c_fit - c_ref, axis=1) / diam
    c_all, _ = centers_rots(
        ref.quaternions_wxyz, ref.translations, np.arange(len(ref_names))
    )
    diam_all = np.max(np.linalg.norm(c_all[:, None, :] - c_all[None, :, :], axis=2))
    f_ref = ref.cameras[0].focal_lengths[0]
    print(
        f"vs reference {ref_files[0].name} ({len(common)} common images; "
        f"subset spans {100 * diam / diam_all:.0f}% of the reference rig):"
    )
    print(
        f"  camera rotation err: mean {rot_err.mean():.2f}, "
        f"median {np.median(rot_err):.2f}, max {rot_err.max():.2f} deg"
    )
    print(
        f"  camera center err:   mean {100 * cen_err.mean():.2f}%, "
        f"median {100 * np.median(cen_err):.2f}%, "
        f"max {100 * cen_err.max():.2f}% of subset diameter"
    )
    print(
        f"  focal: fast {f_est:.1f} px vs reference {f_ref:.1f} px "
        f"({ref.cameras[0].to_dict()['model']}) — "
        f"{100 * (f_est / f_ref - 1):+.1f}%"
    )


# ── Main ─────────────────────────────────────────────────────────────────────


def main():
    global _CAM_WH, _RANK_O
    data = load_clusters()
    obs_c, obs_i, u = data["obs_c"], data["obs_i"], data["obs_uv"]
    n_img, n_cl = data["n_img"], data["n_cl"]
    _CAM_WH = tuple(data["dims"][0])
    _RANK_O = data["adm_rank"][obs_c]
    bam = _RANK_O < MAX_CL
    print(
        f"{WS}: {n_img} images, {n_cl} clusters "
        f"(best {min(MAX_CL, n_cl)} in the BAs), "
        f"{len(obs_c)} observations [{elapsed():.1f}s]"
    )

    # Probe growth at a nominal focal: the core is near-affine, so growth
    # succeeds across a wide f range and its only job here is to build the
    # core geometry and resolve the reflection hypothesis (the two grown
    # candidates' inlier fractions).  The focal itself is chosen later, on
    # geometry wide enough to observe it.
    #
    # Seed groups are tried in covisibility order with a PARALLAX gate: a
    # video's most-mutually-covisible frames are where the camera moved
    # least, and a near-static core fits any focal at high inlier fraction
    # while its depths are unusable (DinoLedge seeded on a static clip).
    f_probe = 0.9 * max(_CAM_WH)
    cap = min(n_img, SCAN_CAP)
    covis = build_covisibility(obs_c, obs_i, n_img, n_cl)
    groups = list(itertools.islice(covis.seed_groups(), 8))
    best = None
    for chunk in range(0, len(groups), 2):
        seed = []
        for group in groups[chunk : chunk + 2]:
            imgs = np.asarray(group)
            wd = factorize_window(obs_c, obs_i, u, imgs)
            seed.append((imgs, wd))
            state = "sparse" if wd is None else f"{len(wd[2][2])} span-2 clusters"
            print(f"seed group {[int(k) for k in imgs]}: {state} [{elapsed():.1f}s]")
        cand = None
        for hyp, grown in enumerate(
            grow_to_cap(seed, f_probe, obs_c, obs_i, u, n_img, n_cl, cap, bam)
        ):
            rvec, tvec, pts, posed, med_inl = grown
            live = ba_rows(posed[obs_i] & ~np.isnan(pts[obs_c, 0]) & bam, obs_i)
            _, rvec, tvec, pts, res, _ = bundle_adjust(
                obs_c[live],
                obs_i[live],
                u[live],
                rvec,
                tvec,
                pts,
                f_probe,
                n_img,
                n_cl,
                opt_f=False,
                schedule=((12.0, 2.0), (4.0, 1.0)),
                max_nfev=30,
            )
            denom = ba_rows(posed[obs_i] & bam, obs_i)
            inl = float((res < 2.0).sum() / max(int(denom.sum()), 1))
            par = core_parallax(rvec, tvec, pts, posed, obs_c, obs_i)
            print(
                f"probe hyp {hyp}: poses {int(posed.sum())}/{n_img}, "
                f"inlier<2px {100 * inl:5.1f}%, parallax {par:.2f} deg "
                f"[{elapsed():.1f}s]"
            )
            if cand is None or inl > cand[0]:
                cand = (inl, par, rvec, tvec, pts, posed, med_inl)
        if cand is None:
            continue
        if best is None or cand[1] > best[1]:
            best = cand
        if cand[1] >= 1.0:  # enough parallax to trust the depths
            break
        print("core parallax too low (near-static seed); trying next groups")
    if best is None:
        raise SystemExit("no seed group produced a reconstruction")

    _, par, rvec, tvec, pts, posed, med_inl = best
    print(f"probe done (parallax {par:.2f} deg) [{elapsed():.1f}s]; widening")
    # Ladder-widen for focal observability: a sliver of a long orbit fits
    # ANY focal at high inlier fraction (bas-relief compensation), so both
    # the focal scan and the release only mean something on a wide arc.
    rvec, tvec, pts, posed = widen(
        rvec, tvec, pts, posed, f_probe, obs_c, obs_i, u, n_img, n_cl, bam,
        gate=0.35 * med_inl,
    )
    print(f"widened to {int(posed.sum())}/{n_img} images [{elapsed():.1f}s]")

    # Photometric verification of the widened set: anchors (valid points,
    # >= 3 posed views, best stored warp-consistency first) are localized
    # in every posed view; an image that keeps too little photometric
    # support is a junk rung the geometric gates missed — un-pose it before
    # the focal scan.  The verified keypoints themselves stay OUT of the
    # estimation: the localization renders through the current geometry,
    # so they carry a bias toward the probe focal (seoul's scan moved from
    # 336 to 432 when fed them) — appearance VERIFIES, geometry decides.
    pv = np.bincount(obs_c[posed[obs_i]], minlength=n_cl)
    cand_a = np.nonzero(~np.isnan(pts[:, 0]) & (pv >= 3))[0]
    cand_a = cand_a[np.argsort(data["cl_quality"][cand_a], kind="stable")]
    anchors = cand_a[:N_ANCHORS]
    a_of_cl = np.full(n_cl, -1, np.int64)
    a_of_cl[anchors] = np.arange(len(anchors))
    sub = np.nonzero(posed)[0]
    rows = np.nonzero((a_of_cl[obs_c] >= 0) & posed[obs_i])[0]
    ai = a_of_cl[obs_c[rows]]
    order = np.argsort(ai, kind="stable")
    rows, ai = rows[order], ai[order]
    sub_of_full = np.full(n_img, -1, np.int64)
    sub_of_full[sub] = np.arange(len(sub))
    a_c, a_i, _a_uv = localize_anchors(
        data["names"],
        sub,
        rvec,
        tvec,
        f_probe,
        pts[anchors],
        ai,
        sub_of_full[obs_i[rows]],
        data["obs_f"][rows],
    )
    kept_per_img = np.bincount(a_i, minlength=n_img)
    floor = max(15, 0.2 * np.median(kept_per_img[sub]))
    for j in sub:
        if kept_per_img[j] < floor:
            posed[j] = False
            print(f"  un-posing image {j}: {kept_per_img[j]} verified obs "
                  f"(floor {floor:.0f})")
    print(f"photometric verify kept {int(posed.sum())}/{len(sub)} images "
          f"[{elapsed():.1f}s]")
    denom = ba_rows(posed[obs_i] & bam, obs_i)

    # Focal scan on the widened geometry: per candidate, rescale the
    # translations (depth scale ~ f), retriangulate, staged fixed-f BA.
    # Two phases: a capped-iteration pass ranks all candidates cheaply, and
    # heavier refits decide between the top two (neighbouring candidates
    # can rank within a point of each other — DinoLedge flips at 52 vs 54 —
    # and the light pass is not reliable at that margin).
    def scan_candidate(f_try, nfev):
        scale = f_try / f_probe
        rv_t, tv_t = rvec.copy(), tvec * scale
        rot = Rotation.from_rotvec(rv_t).as_matrix()
        p_t = triangulate(obs_c, obs_i, u, rot, tv_t, posed, n_cl, f_try)
        live = ba_rows(posed[obs_i] & ~np.isnan(p_t[obs_c, 0]) & bam, obs_i)
        _, rv_t, tv_t, p_t, res, _ = bundle_adjust(
            obs_c[live],
            obs_i[live],
            u[live],
            rv_t,
            tv_t,
            p_t,
            f_try,
            n_img,
            n_cl,
            opt_f=False,
            max_nfev=nfev,
        )
        inl = float((res < 2.0).sum() / max(int(denom.sum()), 1))
        return inl, f_try, rv_t, tv_t, p_t

    coarse = []
    for f_try in np.asarray(F_GRID) * max(_CAM_WH):
        cand = scan_candidate(f_try, 25)
        coarse.append(cand)
        print(f"f={cand[1]:6.1f}: inlier<2px {100 * cand[0]:5.1f}% "
              f"[{elapsed():.1f}s]")
    coarse.sort(key=lambda t: -t[0])
    finals = [scan_candidate(c[1], 60) for c in coarse[:2]]
    for c in finals:
        print(f"f={c[1]:6.1f} (refit): inlier<2px {100 * c[0]:5.1f}% "
              f"[{elapsed():.1f}s]")
    best = max(finals, key=lambda t: t[0])

    inl0, f, rvec, tvec, pts = best
    print(f"scan winner: f = {f:.1f} [{elapsed():.1f}s]; releasing f")
    # Iterated release: full schedule (the wide first trim + inter-round
    # retriangulation is what lets f keep walking — the structure absorbs a
    # wrong f and must be re-formed as f moves).  Stop when f stabilizes;
    # keep the best-fit state seen.
    inl, f_prev = inl0, f
    kept = (inl0, f, rvec, tvec, pts)
    for _ in range(3):
        live = ba_rows(posed[obs_i] & ~np.isnan(pts[obs_c, 0]) & bam, obs_i)
        f, rvec, tvec, pts, res, _ = bundle_adjust(
            obs_c[live],
            obs_i[live],
            u[live],
            rvec,
            tvec,
            pts,
            f,
            n_img,
            n_cl,
            opt_f=True,
            max_nfev=30,
        )
        inl = float((res < 2.0).sum() / max(int(denom.sum()), 1))
        # The affine collapse is an UPWARD escape: on narrow or shallow
        # geometry, f -> inf keeps fitting better (rising inlier fraction —
        # which also means the keep-best rule cannot be trusted upward), so
        # a release drifting above the scan winner contradicts the scan's
        # trim-consistent ranking and loses to it (every legitimate upward
        # walk observed stayed under +10%).  Downward walks are legitimate
        # (a true focal below the grid floor releases past the smallest
        # candidate) up to a plausibility floor.
        if f > 1.15 * best[1] or f < 0.3 * max(_CAM_WH):
            print(f"release left the scan basin (f = {f:.0f}); keeping previous")
            break
        if inl > kept[0]:
            kept = (inl, f, rvec, tvec, pts)
        if abs(f - f_prev) < 0.01 * f_prev:
            break
        f_prev = f
    inl, f, rvec, tvec, pts = kept
    print(
        f"\nfast pinhole estimate: f = {f:.1f} px on {int(posed.sum())}/{n_img} "
        f"images, inlier<2px {100 * inl:.1f}% [{elapsed():.1f}s]"
    )
    compare_to_reference(data["names"], rvec, tvec, f, posed)

    out = {
        "focal_px": float(f),
        "width": _CAM_WH[0],
        "height": _CAM_WH[1],
        "posed_images": [n for j, n in enumerate(data["names"]) if posed[j]],
        "rvec": rvec[posed].tolist(),
        "tvec": tvec[posed].tolist(),
        "elapsed_s": round(elapsed(), 2),
    }
    (WS / "fast-pinhole.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {WS / 'fast-pinhole.json'}")


if __name__ == "__main__":
    main()
