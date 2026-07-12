# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment: pinhole-only coarse reconstruction from cluster patches.

Starting from a workspace holding images and a `*-clusters-patches.matches`
file (sift extraction -> cluster matching -> cluster-patches), and using only
a pinhole camera model, bootstrap a coarse 3D reconstruction and write it to
a `.sfmr` file — no COLMAP solver involved.

Pipeline:
  1. Load patch clusters; refined member positions are read directly from
     the stored affines' last column (`member_affines[k][:, 2]` holds the
     absolute keypoint position since .matches format version 4), and the
     image dimensions from the images section — no per-image .sift reads.
  2. Group images by cluster covisibility (shared-cluster counts) — no
     sequence order is assumed.  Affine (weak-perspective) ALS factorization
     of candidate seed groups (a single global factorization breaks on wide
     baselines) + Tomasi–Kanade metric upgrade, both reflection hypotheses.
  3. Seed a perspective solve on the best group (a small fixed-focal BA
     also resolves the reflection), then grow incrementally: the
     next-best-view image (most observations of valid points) is resected
     pose-only against the global structure (trimmed iterations, most-
     covisible posed poses as inits), new clusters are triangulated as
     they gain posed views, short global BAs run every few images.
  4. Steps 2–3 run per candidate focal on a small grid with f held FIXED —
     the focal is unobservable from a weak init (the residual decreases
     monotonically toward the affine limit), but with a converged geometry
     the inlier fraction peaks near the true focal.  The scan caps growth
     at ~20 images; the winner grows fully and its BA then releases f.
  5. Report reprojection stats and, when a reference solve exists in the
     workspace, camera errors after similarity alignment; save the result
     as `sfmr/bootstrap-pinhole.sfmr`.

Run: pixi run -e dev python scripts/exp_pinhole_bootstrap.py <workspace> [ref.sfmr]

The optional second argument names the reference solve to compare against
(it may live in another workspace, e.g. a full-sequence solve when
bootstrapping a frame subset — images are matched by workspace-relative
name).  Default: the first non-bootstrap .sfmr in the workspace.
"""

import sys
from pathlib import Path

import numpy as np
import scipy.optimize
import scipy.sparse
from scipy.spatial.transform import Rotation

WS = Path(sys.argv[1] if len(sys.argv) > 1 else "e_seoul_ws")
REF = Path(sys.argv[2]) if len(sys.argv) > 2 else None
MIN_SPAN_BA = 2  # min distinct images for a cluster to become a point
MAX_CLUSTERS = 10000  # cap for the scipy BAs; highest-span clusters win
F_GRID = [0.55, 0.7, 0.9, 1.2, 1.6]  # focal candidates, in units of max(w, h)
TRIM_PX = 4.0  # BA inter-round observation trim threshold


# ── Data loading ─────────────────────────────────────────────────────────────


def load_clusters():
    """Patch clusters as flat observation arrays with refined positions.

    Everything geometric comes straight from the .matches file: image
    dimensions from the images section and member positions from the stored
    affines' last column (the absolute refined keypoint position).
    """
    from sfmtool._sfmtool.io import read_matches

    patches = sorted(WS.glob("matches/*-clusters-patches.matches"))
    data = read_matches(patches[0])
    names = list(data["image_names"])
    dims = [(int(w), int(h)) for w, h in np.asarray(data["image_dims"])]

    starts = np.asarray(data["cluster_starts"])
    mi = np.asarray(data["member_images"])
    mf = np.asarray(data["member_features"])
    st = np.asarray(data["member_status"])
    refs = np.asarray(data["reference_members"])
    aff = np.asarray(data["member_affines"])

    # First pass: member selections and spans of every usable cluster.
    usable = []
    for c in range(len(starts) - 1):
        lo, hi = int(starts[c]), int(starts[c + 1])
        if refs[c] == np.iinfo(np.uint32).max:
            continue
        sel = np.nonzero((st[lo:hi] == 0) | (st[lo:hi] == 1))[0] + lo
        span = len(np.unique(mi[sel]))
        if span >= MIN_SPAN_BA:
            usable.append((span, c, sel))

    # Cap for the scipy BAs: keep the highest-span clusters (deterministic;
    # ties broken by cluster id).  A coarse bootstrap doesn't need them all.
    if len(usable) > MAX_CLUSTERS:
        usable.sort(key=lambda t: (-t[0], t[1]))
        dropped = len(usable) - MAX_CLUSTERS
        usable = sorted(usable[:MAX_CLUSTERS], key=lambda t: t[1])
        print(f"capped clusters: kept {MAX_CLUSTERS} by span, dropped {dropped}")

    obs_c, obs_i, obs_f, obs_uv = [], [], [], []
    n_cl = 0
    for _span, c, sel in usable:
        for k in sel:
            obs_c.append(n_cl)
            obs_i.append(int(mi[k]))
            obs_f.append(int(mf[k]))
            # The affine's last column is the member's absolute refined
            # keypoint position (identity | x_ref for the reference row).
            obs_uv.append(aff[k, :, 2])
        n_cl += 1

    return {
        "names": names,
        "dims": dims,
        "obs_c": np.asarray(obs_c),
        "obs_i": np.asarray(obs_i),
        "obs_f": np.asarray(obs_f),
        "obs_uv": np.asarray(obs_uv, dtype=np.float64),
        "n_img": len(names),
        "n_cl": n_cl,
    }


# ── Affine factorization (ALS with trimming) ─────────────────────────────────


def als_factorize(obs_c, obs_i, u, n_img, n_cl, rounds=25, trim=0.05):
    grid = np.zeros((2 * n_img, n_cl))
    cnt = np.zeros((n_img, n_cl))
    for c, i, uv in zip(obs_c, obs_i, u):
        grid[2 * i : 2 * i + 2, c] = uv
        cnt[i, c] = 1
    row_mean = grid.sum(axis=1, keepdims=True) / np.maximum(
        np.repeat(cnt.sum(axis=1), 2)[:, None], 1
    )
    filled = np.where(np.repeat(cnt, 2, axis=0) > 0, grid, row_mean)
    filled -= row_mean
    _, _, vt = np.linalg.svd(filled, full_matrices=False)
    x_pts = vt[:3].T

    keep = np.ones(len(obs_c), bool)
    m_cam = np.zeros((n_img, 2, 3))
    t_cam = np.zeros((n_img, 2))
    for it in range(rounds):
        for i in range(n_img):
            s = keep & (obs_i == i)
            if s.sum() < 4:
                continue
            xh = np.concatenate([x_pts[obs_c[s]], np.ones((s.sum(), 1))], axis=1)
            sol = np.linalg.lstsq(xh, u[s], rcond=None)[0]
            m_cam[i] = sol[:3].T
            t_cam[i] = sol[3]
        for c in range(n_cl):
            s = keep & (obs_c == c)
            if s.sum() < 2:
                continue
            a = m_cam[obs_i[s]].reshape(-1, 3)
            b = (u[s] - t_cam[obs_i[s]]).reshape(-1)
            x_pts[c] = np.linalg.lstsq(a, b, rcond=None)[0]
        res = u - (np.einsum("nij,nj->ni", m_cam[obs_i], x_pts[obs_c]) + t_cam[obs_i])
        if trim > 0 and it >= rounds // 2:
            rn = np.linalg.norm(res, axis=1)
            keep = rn < np.quantile(rn[keep], 1 - trim)
    return m_cam, t_cam, x_pts, res, keep


# ── Metric upgrade (Tomasi–Kanade) ───────────────────────────────────────────


def metric_upgrade(m_cam, used):
    """Solve the 3x3 gauge A so that rows of M_i·A are orthogonal and equal
    norm.  Linear least squares on the symmetric Q = A·Aᵀ (6 unknowns).
    Returns both reflection hypotheses of A."""
    rows_a, rows_b = [], []

    def sym_row(p, q):
        # coefficients of pᵀ Q q over Q's 6 upper-triangle entries
        return np.array(
            [
                p[0] * q[0],
                p[0] * q[1] + p[1] * q[0],
                p[0] * q[2] + p[2] * q[0],
                p[1] * q[1],
                p[1] * q[2] + p[2] * q[1],
                p[2] * q[2],
            ]
        )

    for i in np.nonzero(used)[0]:
        m1, m2 = m_cam[i]
        rows_a.append(sym_row(m1, m1) - sym_row(m2, m2))
        rows_a.append(sym_row(m1, m2))
        rows_b.append(sym_row(m1, m1) + sym_row(m2, m2))
    a = np.asarray(rows_a)
    # Normalization: mean row-norm² = 1 (avoids the trivial Q = 0).
    a = np.vstack([a, np.mean(rows_b, axis=0)[None, :]])
    b = np.zeros(len(a))
    b[-1] = 2.0
    qv = np.linalg.lstsq(a, b, rcond=None)[0]
    q_mat = np.array(
        [
            [qv[0], qv[1], qv[2]],
            [qv[1], qv[3], qv[4]],
            [qv[2], qv[4], qv[5]],
        ]
    )
    w, v = np.linalg.eigh(q_mat)
    w = np.maximum(w, 1e-8 * w.max())
    a_up = v @ np.diag(np.sqrt(w))
    return a_up, a_up @ np.diag([1.0, 1.0, -1.0])


def weak_perspective_poses(m_cam, t_cam, a_up, used):
    """Per-image rotation + scale from the metric-upgraded affine cameras."""
    n_img = len(m_cam)
    rot = np.zeros((n_img, 3, 3))
    scale = np.zeros(n_img)
    for i in range(n_img):
        if not used[i]:
            continue
        m = m_cam[i] @ a_up
        s = 0.5 * (np.linalg.norm(m[0]) + np.linalg.norm(m[1]))
        r3 = np.cross(m[0], m[1])
        stack = np.vstack([m / s, r3 / max(np.linalg.norm(r3), 1e-12)])
        uu, _, vv = np.linalg.svd(stack)
        r = uu @ vv
        if np.linalg.det(r) < 0:
            r = uu @ np.diag([1, 1, -1]) @ vv
        rot[i] = r
        scale[i] = s
    return rot, scale


# ── Covisibility grouping ────────────────────────────────────────────────────
#
# No sequence order is assumed: the natural grouping is how many clusters a
# pair of images shares.  High mutual covisibility implies nearby viewpoints,
# which is exactly what the weak-perspective factorization needs from a seed
# group, and the same counts drive the growth order and the resection inits.


def covisibility(obs_c, obs_i, n_img):
    """W[i, j] = number of clusters observed in both image i and image j."""
    w = np.zeros((n_img, n_img), dtype=np.int32)
    order = np.argsort(obs_c, kind="stable")
    oc, oi = obs_c[order], obs_i[order]
    for imgs in np.split(oi, np.nonzero(np.diff(oc))[0] + 1):
        uu = np.unique(imgs)
        w[np.ix_(uu, uu)] += 1
    np.fill_diagonal(w, 0)
    return w


def pick_seed_groups(w, size=5, count=2, min_shared=8):
    """Up to ``count`` disjoint candidate seed groups: greedily grow from a
    strongest remaining edge, each step adding the image with the best
    minimum shared-cluster count against the whole group (mutual
    covisibility, not a hub-and-spokes star)."""
    w = w.copy()
    groups = []
    for _ in range(count):
        if w.max() < min_shared:
            break
        i, j = np.unravel_index(np.argmax(w), w.shape)
        group = [int(i), int(j)]
        while len(group) < size:
            mins = w[:, group].min(axis=1)
            mins[group] = -1
            k = int(np.argmax(mins))
            if mins[k] < min_shared:
                break
            group.append(k)
        groups.append(np.array(sorted(group)))
        w[group, :] = 0
        w[:, group] = 0
    return groups


def window_spans(obs_c, obs_i, u, imgs, min_span):
    """Observation selection for clusters seen in >= min_span window images."""
    inw = np.isin(obs_i, imgs)
    cl, il = obs_c[inw], np.searchsorted(imgs, obs_i[inw])
    span = np.zeros(cl.max() + 1 if len(cl) else 1, int)
    for c in np.unique(cl):
        span[c] = len(np.unique(il[cl == c]))
    good_cl = np.nonzero(span >= min_span)[0]
    sel = inw.copy()
    sel[inw] = np.isin(cl, good_cl)
    uniq, c2 = np.unique(obs_c[sel], return_inverse=True)
    return sel, np.searchsorted(imgs, obs_i[sel]), uniq, c2


def factorize_window(obs_c, obs_i, u, imgs, min_span=3):
    """Factorize the clusters seen in >= min_span images of the window.

    Returns (both metric hypotheses as (rot, scale, t_aff) in the window's
    local frame, used mask, span-2 selection for the window mini-BA) or None
    when the window is too sparse.
    """
    sel, il, uniq, c2 = window_spans(obs_c, obs_i, u, imgs, min_span)
    if sel.sum() < 30:
        return None
    m_cam, t_aff, _, _, keep = als_factorize(c2, il, u[sel], len(imgs), len(uniq))
    used = np.bincount(il[keep], minlength=len(imgs)) >= 4
    hyps = []
    for a_up in metric_upgrade(m_cam, used):
        rot, scale = weak_perspective_poses(m_cam, t_aff, a_up, used)
        if (scale[used] > 0).all():
            hyps.append((rot, scale, t_aff))
    ba_sel = window_spans(obs_c, obs_i, u, imgs, 2)
    return hyps, used, ba_sel


def pose_refine(uv, x_pts, rv0, tv0, f):
    """Pose-only resection of one image against known 3D points.

    Trimmed iterations: repeatedly refit L2 on the best-fitting 60% of the
    observations.  A plain L2 warm-up is dragged by the junk observations'
    leverage, and a robust loss has near-zero gradient when every residual
    starts as a 100 px "outlier" — trimming from a decent init has neither
    problem.  A final refit uses the < 3 px inliers."""

    def rn_of(x):
        xc = Rotation.from_rotvec(x[:3]).apply(x_pts) + x[3:]
        z = np.maximum(xc[:, 2], 1e-6)
        return np.linalg.norm(f * xc[:, :2] / z[:, None] - uv, axis=1)

    def fit(x0, mask):
        pts_sel, uv_sel = x_pts[mask], uv[mask]

        def resid(x):
            xc = Rotation.from_rotvec(x[:3]).apply(pts_sel) + x[3:]
            z = np.maximum(xc[:, 2], 1e-6)
            return (f * xc[:, :2] / z[:, None] - uv_sel).ravel()

        return scipy.optimize.least_squares(resid, x0, max_nfev=60).x

    x0 = np.concatenate([rv0, tv0])
    rn = rn_of(x0)
    for _ in range(5):
        x0 = fit(x0, rn <= np.quantile(rn, 0.6))
        rn = rn_of(x0)
    inl = rn < 3.0
    if inl.sum() >= 6:
        x0 = fit(x0, inl)
        rn = rn_of(x0)
    return x0[:3], x0[3:], float((rn < 3.0).mean())


def fill_new_points(pts, obs_c, obs_i, u, rvec, tvec, posed, f):
    """DLT-triangulate clusters that lack a point but now have >= 2 posed
    observations.  Existing points are left untouched."""
    need = np.isnan(pts[:, 0])[obs_c] & posed[obs_i]
    if not need.any():
        return pts
    uniq, c2 = np.unique(obs_c[need], return_inverse=True)
    rot = Rotation.from_rotvec(rvec).as_matrix()
    newp = triangulate(c2, obs_i[need], u[need], rot, tvec, posed, len(uniq), f)
    out = pts.copy()
    out[uniq] = newp
    return out


def grow_reconstruction(
    grp_data, f0, obs_c, obs_i, u, n_img, n_cl, covis, max_images=None
):
    """Incremental bootstrap, no sequence order assumed.

    Seed a perspective solve on a candidate covisibility group (all groups x
    both reflection hypotheses, best inlier fraction wins), then repeatedly
    pose the unposed image with the most observations of valid 3D points
    (next-best-view) by trimmed pose-only resection — initialised from its
    most-covisible posed images — triangulating new clusters as they gain
    posed views, with a short global BA every few images.

    ``max_images`` caps growth (used by the focal scan, which only needs
    enough geometry to rank candidates).  Returns (rvec, tvec, pts, posed)
    or None.
    """
    grow_schedule = [(30.0, 3.0), (8.0, 1.5)]
    ba_every = max(3, min(8, n_img // 10))

    best = None
    for imgs, wd in grp_data:
        if wd is None:
            continue
        hyps, used, (sel, il, cl_ids, c2) = wd
        for rot0, scale, t_aff in hyps:
            trans0 = perspective_init(rot0, scale, t_aff, used, f0)
            pts_w = triangulate(c2, il, u[sel], rot0, trans0, used, len(cl_ids), f0)
            ok = ~np.isnan(pts_w[:, 0])[c2] & used[il]
            _, rvw, tvw, p_w, _, _, inl = bundle_adjust(
                c2[ok],
                il[ok],
                u[sel][ok],
                rot0,
                trans0,
                pts_w,
                f0,
                len(imgs),
                len(cl_ids),
                opt_f=False,
                verbose=False,
            )
            if best is None or inl > best[0]:
                best = (inl, imgs, used, cl_ids, rvw, tvw, p_w)
    if best is None:
        return None
    _, imgs, used, cl_ids, rvw, tvw, p_w = best

    rvec = np.zeros((n_img, 3))
    tvec = np.tile([0.0, 0.0, f0], (n_img, 1))
    posed = np.zeros(n_img, bool)
    pts = np.full((n_cl, 3), np.nan)
    for k, i in enumerate(imgs):
        if used[k]:
            rvec[i], tvec[i], posed[i] = rvw[k], tvw[k], True
    pts[cl_ids] = p_w

    since_ba = 0
    while max_images is None or posed.sum() < max_images:
        # Next-best-view: most observations of currently-valid points.
        cand = ~posed[obs_i] & ~np.isnan(pts[obs_c, 0])
        if not cand.any():
            break
        cnt = np.bincount(obs_i[cand], minlength=n_img)
        i = int(np.argmax(cnt))
        if cnt[i] < 6:
            break
        s = (obs_i == i) & ~np.isnan(pts[obs_c, 0])
        # Resect from the most-covisible posed images' poses; keep the best.
        posed_idx = np.nonzero(posed)[0]
        inits = posed_idx[np.argsort(-covis[i, posed_idx])][:3]
        found = None
        for j in inits:
            rv, tv, inl = pose_refine(u[s], pts[obs_c[s]], rvec[j], tvec[j], f0)
            if found is None or inl > found[0]:
                found = (inl, rv, tv)
            if inl > 0.4:
                break
        _, rvec[i], tvec[i] = found
        posed[i] = True
        pts = fill_new_points(pts, obs_c, obs_i, u, rvec, tvec, posed, f0)
        since_ba += 1
        if since_ba >= ba_every:
            since_ba = 0
            live = posed[obs_i] & ~np.isnan(pts[obs_c, 0])
            rot = Rotation.from_rotvec(rvec).as_matrix()
            ba = bundle_adjust(
                obs_c[live],
                obs_i[live],
                u[live],
                rot,
                tvec,
                pts,
                f0,
                n_img,
                n_cl,
                opt_f=False,
                verbose=False,
                schedule=grow_schedule,
            )
            rvec, tvec, pts = ba[1], ba[2], ba[3]
    return rvec, tvec, pts, posed


# ── Perspective conversion + triangulation ───────────────────────────────────


def perspective_init(rot, scale, t_cam, used, f0):
    """Weak-perspective -> pinhole poses: depth of the point centroid is
    f0/s, lateral offset t/s (COLMAP-style camera frame, +z forward)."""
    trans = np.zeros((len(rot), 3))
    for i in np.nonzero(used)[0]:
        trans[i] = [t_cam[i, 0] / scale[i], t_cam[i, 1] / scale[i], f0 / scale[i]]
    return trans


def triangulate(obs_c, obs_i, u, rot, trans, used, n_cl, f):
    """Linear DLT triangulation of every cluster from the posed images."""
    pts = np.full((n_cl, 3), np.nan)
    ok_img = used[obs_i]
    for c in range(n_cl):
        s = (obs_c == c) & ok_img
        if s.sum() < 2:
            continue
        rows = []
        for i, uv in zip(obs_i[s], u[s]):
            p = np.hstack([rot[i], trans[i][:, None]])  # 3x4, cam frame
            rows.append(uv[0] * p[2] - f * p[0])
            rows.append(uv[1] * p[2] - f * p[1])
        m = np.asarray(rows)
        _, _, vt = np.linalg.svd(m)
        h = vt[-1]
        if abs(h[3]) < 1e-12:
            continue
        pts[c] = h[:3] / h[3]
    return pts


# ── Bundle adjustment ────────────────────────────────────────────────────────


def project(f, rvec, tvec, pts, obs_i, obs_c):
    xc = Rotation.from_rotvec(rvec[obs_i]).apply(pts[obs_c]) + tvec[obs_i]
    z = np.maximum(xc[:, 2], 1e-6)
    return f * xc[:, :2] / z[:, None], xc[:, 2]


def solve_round(obs_c, obs_i, u, f, rvec, tvec, pts, n_img, f_scale, opt_f):
    """One robust sparse least-squares pass over the given observations.

    Images and points are compact-indexed over what the observations touch;
    untouched poses/points pass through unchanged."""
    uniq, oc2 = np.unique(obs_c, return_inverse=True)
    uimg, oi2 = np.unique(obs_i, return_inverse=True)
    p0 = pts[uniq]
    n_pt, n_im = len(uniq), len(uimg)
    nf = 1 if opt_f else 0

    def unpack(x):
        fv = x[0] if opt_f else f
        rv = x[nf : nf + 3 * n_im].reshape(-1, 3)
        tv = x[nf + 3 * n_im : nf + 6 * n_im].reshape(-1, 3)
        pv = x[nf + 6 * n_im :].reshape(-1, 3)
        return fv, rv, tv, pv

    def resid(x):
        fv, rv, tv, pv = unpack(x)
        proj, _ = project(fv, rv, tv, pv, oi2, oc2)
        return (proj - u).ravel()

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
        [[f] if opt_f else [], rvec[uimg].ravel(), tvec[uimg].ravel(), p0.ravel()]
    )
    sol = scipy.optimize.least_squares(
        resid,
        x0,
        jac_sparsity=spar,
        loss="soft_l1",
        f_scale=f_scale,
        max_nfev=60,
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
    rot,
    trans,
    pts,
    f0,
    n_img,
    n_cl,
    opt_f,
    verbose=True,
    schedule=None,
):
    """Staged robust BA: trim gross outliers and behind-camera observations
    before each solve, re-triangulate every cluster from the refined cameras
    between rounds (re-admitting points the bad init lost)."""
    rvec = Rotation.from_matrix(rot).as_rotvec()
    tvec = trans.copy()
    f = f0
    all_img = np.ones(n_img, bool)
    if schedule is None:
        schedule = [(50.0, 5.0), (12.0, 2.0), (TRIM_PX, 1.0)]

    for rnd, (thresh, f_scale) in enumerate(schedule):
        if rnd > 0:
            rot_now = Rotation.from_rotvec(rvec).as_matrix()
            pts = triangulate(obs_c, obs_i, u, rot_now, tvec, all_img, n_cl, f)
        proj, depth = project(f, rvec, tvec, pts, obs_i, obs_c)
        rn = np.linalg.norm(proj - u, axis=1)
        with np.errstate(invalid="ignore"):
            keep = (rn < thresh) & (depth > 1e-3 * f) & ~np.isnan(pts[obs_c, 0])
        surv = np.bincount(obs_c[keep], minlength=n_cl)
        keep &= surv[obs_c] >= MIN_SPAN_BA
        f, rvec, tvec, pts = solve_round(
            obs_c[keep],
            obs_i[keep],
            u[keep],
            f,
            rvec,
            tvec,
            pts,
            n_img,
            f_scale,
            opt_f,
        )
        proj, depth = project(f, rvec, tvec, pts, obs_i, obs_c)
        rn = np.linalg.norm(proj - u, axis=1)
        if verbose:
            print(
                f"  BA round {rnd} (trim {thresh:.0f}px): f {f:.1f}, "
                f"median reproj {np.nanmedian(rn):.2f} px on {keep.sum()} obs"
            )

    with np.errstate(invalid="ignore"):
        keep = (rn < TRIM_PX) & (depth > 1e-3 * f) & ~np.isnan(pts[obs_c, 0])
    surv = np.bincount(obs_c[keep], minlength=n_cl)
    keep &= surv[obs_c] >= MIN_SPAN_BA
    res = np.where(np.isnan(rn), np.inf, rn)
    inlier2 = float((res < 2.0).mean())
    return f, rvec, tvec, pts, keep, res, inlier2


# ── Evaluation against a reference solve ─────────────────────────────────────


def umeyama(src, dst):
    """Similarity (s, R, t) minimizing ||s·R·src + t − dst||, no reflection."""
    mu_s, mu_d = src.mean(axis=0), dst.mean(axis=0)
    sc, dc = src - mu_s, dst - mu_d
    cov = dc.T @ sc / len(src)
    uu, dd, vv = np.linalg.svd(cov)
    sgn = np.sign(np.linalg.det(uu @ vv))
    r = uu @ np.diag([1, 1, sgn]) @ vv
    var = (sc**2).sum() / len(src)
    s = (dd * [1, 1, sgn]).sum() / var
    t = mu_d - s * r @ mu_s
    return s, r, t


def compare_to_reference(names, rvec, tvec, f_est):
    """Compare against the first non-bootstrap solve in the workspace.

    Our BA poses are COLMAP-convention; the reference ``.sfmr`` is canonical.
    Convert ours to canonical first so the per-camera rotation errors are
    meaningful (the world-frame difference is absorbed by the alignment).
    """
    if REF is not None:
        ref_files = [REF]
    else:
        ref_files = sorted(
            p for p in WS.glob("sfmr/*.sfmr") if p.name != "bootstrap-pinhole.sfmr"
        )
    if not ref_files:
        print("no reference solve found; skipping comparison")
        return
    from sfmtool._sfmtool import SfmrReconstruction
    from sfmtool.colmap.io import _colmap_poses_points_to_canonical

    q_colmap = Rotation.from_rotvec(rvec).as_quat()[:, [3, 0, 1, 2]]
    q_wxyz, t_xyz, _ = _colmap_poses_points_to_canonical(
        q_colmap, tvec, np.zeros((0, 3)), apply_world_rotation=True
    )

    ref = SfmrReconstruction.load(ref_files[0])
    ref_names = list(ref.image_names)
    common = [n for n in names if n in ref_names]
    if len(common) < 3:
        print(f"only {len(common)} common images with {ref_files[0].name}; skipping")
        return

    def centers_rots(qs, ts, order):
        rs = Rotation.from_quat(np.asarray(qs)[order][:, [1, 2, 3, 0]]).as_matrix()
        cs = -np.einsum("nij,ni->nj", rs, np.asarray(ts)[order])
        return cs, rs

    ei = np.array([names.index(n) for n in common])
    ri = np.array([ref_names.index(n) for n in common])
    c_est, r_est = centers_rots(q_wxyz, t_xyz, ei)
    c_ref, r_ref = centers_rots(ref.quaternions_wxyz, ref.translations, ri)

    s, r_al, t_al = umeyama(c_est, c_ref)
    c_fit = (s * (r_al @ c_est.T)).T + t_al
    diam = np.max(np.linalg.norm(c_ref[:, None, :] - c_ref[None, :, :], axis=2))
    cen_err = np.linalg.norm(c_fit - c_ref, axis=1) / diam
    rot_err = Rotation.from_matrix(
        np.einsum("nij,nkj->nik", r_ref, np.einsum("nij,kj->nik", r_est, r_al))
    ).magnitude() * (180 / np.pi)

    cam0 = ref.cameras[0].to_dict()
    f_ref = cam0["parameters"].get("focal_length", cam0["parameters"].get("fx"))
    print(f"\nvs reference {ref_files[0].name} ({len(common)} common images):")
    print(
        f"  camera rotation err: mean {rot_err.mean():.2f}, "
        f"median {np.median(rot_err):.2f}, max {rot_err.max():.2f} deg; "
        f"{(rot_err > 10).sum()} cams > 10 deg"
    )
    print(
        f"  camera center err:   mean {100 * cen_err.mean():.2f}%, "
        f"median {100 * np.median(cen_err):.2f}%, "
        f"max {100 * cen_err.max():.2f}% of scene diameter"
    )
    print(
        f"  focal: bootstrap {f_est:.1f} px vs reference {f_ref:.1f} px "
        f"({cam0['model']})"
    )


# ── Save as .sfmr ────────────────────────────────────────────────────────────


def save_sfmr(data, f, rvec, tvec, pts, keep, res, out_path):
    from sfmtool._sfmtool import SfmrReconstruction
    from sfmtool._sfmtool.geometry import CameraIntrinsics
    from sfmtool._workspace import load_workspace_config
    from sfmtool.colmap.io import (
        _build_sfmr_data_dict,
        _colmap_poses_points_to_canonical,
        _resolve_workspace_and_sift,
        build_metadata,
        finite_positions_xyzw,
    )

    names, dims = data["names"], data["dims"]
    obs_c, obs_i, obs_f = data["obs_c"], data["obs_i"], data["obs_f"]
    w, h = dims[0]

    # Surviving points, renumbered densely; observations grouped by point.
    alive = np.nonzero(np.bincount(obs_c[keep], minlength=len(pts)) >= 2)[0]
    remap = {int(c): k for k, c in enumerate(alive)}
    order = np.argsort(obs_c[keep], kind="stable")
    ko = np.nonzero(keep)[0][order]
    ko = ko[np.isin(obs_c[ko], alive)]

    track_img = obs_i[ko]
    track_feat = obs_f[ko]
    point_idx = np.array([remap[int(c)] for c in obs_c[ko]])
    obs_counts = np.bincount(point_idx, minlength=len(alive))

    positions = pts[alive]
    per_point_err = np.zeros(len(alive), dtype=np.float32)
    np.add.at(per_point_err, point_idx, res[ko].astype(np.float32))
    per_point_err /= np.maximum(obs_counts, 1)

    q_colmap = Rotation.from_rotvec(rvec).as_quat()[:, [3, 0, 1, 2]]
    q_can, t_can, p_can = _colmap_poses_points_to_canonical(
        q_colmap, tvec, positions, apply_world_rotation=True
    )

    (
        workspace_dir,
        _contents,
        resolved_names,
        ft_hashes,
        sc_hashes,
        thumbnails,
    ) = _resolve_workspace_and_sift(names, WS)

    # Colors from the .sift thumbnails at the (scaled) observation position.
    colors = np.zeros((len(alive), 3), dtype=np.uint8)
    uv = data["obs_uv"][ko]
    for k in range(len(ko)):
        th = np.asarray(thumbnails[track_img[k]])
        ty = int(np.clip(uv[k, 1] * th.shape[0] / h, 0, th.shape[0] - 1))
        tx = int(np.clip(uv[k, 0] * th.shape[1] / w, 0, th.shape[1] - 1))
        colors[point_idx[k]] = th[ty, tx]

    camera = CameraIntrinsics.from_dict(
        {
            "model": "SIMPLE_PINHOLE",
            "width": w,
            "height": h,
            "parameters": {
                "focal_length": float(f),
                "principal_point_x": w / 2,
                "principal_point_y": h / 2,
            },
        }
    )

    metadata = build_metadata(
        workspace_dir=workspace_dir,
        output_path=out_path,
        workspace_config=load_workspace_config(workspace_dir),
        operation="cluster_bootstrap",
        tool_name="sfmtool",
        tool_options={"camera_model": "SIMPLE_PINHOLE", "focal_grid": F_GRID},
        image_count=len(names),
        point_count=len(alive),
        observation_count=int(obs_counts.sum()),
        camera_count=1,
    )

    sfmr_dict = _build_sfmr_data_dict(
        cameras=[camera],
        image_names=resolved_names,
        camera_indexes=np.zeros(len(names), dtype=np.uint32),
        quaternions_wxyz=q_can,
        translations_xyz=t_can,
        positions_xyzw=finite_positions_xyzw(p_can),
        colors_rgb=colors,
        reprojection_errors=per_point_err,
        track_image_indexes=track_img,
        track_feature_indexes=track_feat,
        point_indexes=point_idx,
        observation_counts=obs_counts,
        feature_tool_hashes=ft_hashes,
        sift_content_hashes=sc_hashes,
        thumbnails=thumbnails,
        metadata=metadata,
    )

    recon = SfmrReconstruction.from_data(workspace_dir, sfmr_dict)
    recon.save(out_path)
    print(f"\nwrote {out_path} ({len(alive)} points, {int(obs_counts.sum())} obs)")
    return recon


# ── Main ─────────────────────────────────────────────────────────────────────


def main():
    data = load_clusters()
    obs_c, obs_i, u_px = data["obs_c"], data["obs_i"], data["obs_uv"]
    n_img, n_cl = data["n_img"], data["n_cl"]
    dims = np.asarray(data["dims"], dtype=np.float64)
    half = dims[data["obs_i"]] / 2
    u = u_px - half
    print(
        f"{WS}: {n_img} images, {n_cl} clusters (span >= {MIN_SPAN_BA}), "
        f"{len(obs_c)} observations"
    )

    # Covisibility grouping — no sequence order assumed anywhere.
    covis = covisibility(obs_c, obs_i, n_img)
    groups = pick_seed_groups(covis)
    grp_data = []
    for imgs in groups:
        wd = factorize_window(obs_c, obs_i, u, imgs)
        grp_data.append((imgs, wd))
        state = "sparse" if wd is None else f"{len(wd[2][2])} span-2 clusters"
        print(f"seed group {[int(k) for k in imgs]}: {state}")

    # The focal is unobservable at init quality (the residual decreases
    # monotonically toward the affine limit), so grow the reconstruction at
    # each candidate focal with f held FIXED and let the inlier fraction
    # pick.  The seed group's inlier fraction resolves the reflection
    # hypothesis.  The scan caps growth at ~20 images (enough to rank);
    # only the winner grows fully, then releases f.
    f_grid = np.array(F_GRID) * dims.max()
    scan_cap = min(n_img, 20)
    best = None
    for f_try in f_grid:
        grown = grow_reconstruction(
            grp_data,
            f_try,
            obs_c,
            obs_i,
            u,
            n_img,
            n_cl,
            covis,
            max_images=scan_cap,
        )
        if grown is None:
            continue
        g_rvec, g_tvec, pts, posed = grown
        rot = Rotation.from_rotvec(g_rvec).as_matrix()
        ok = posed[obs_i] & ~np.isnan(pts[:, 0])[obs_c]
        ba = bundle_adjust(
            obs_c[ok],
            obs_i[ok],
            u[ok],
            rot,
            g_tvec,
            pts,
            f_try,
            n_img,
            n_cl,
            opt_f=False,
            verbose=False,
        )
        res = ba[5]
        # Rank on ALL observations of the posed subset: clusters that failed
        # to triangulate under this candidate count as misses.
        inl_scan = float((res < 2.0).sum() / max(posed[obs_i].sum(), 1))
        print(
            f"f={f_try:6.1f}: poses {posed.sum()}/{n_img}, "
            f"inlier<2px {100 * inl_scan:5.1f}%, "
            f"median {np.median(res[np.isfinite(res)]):6.2f} px"
        )
        if best is None or inl_scan > best[0]:
            best = (inl_scan, f_try)

    _, f_try = best
    print(f"\nwinner: f = {f_try:.1f}; growing fully, then releasing f")
    grown = grow_reconstruction(grp_data, f_try, obs_c, obs_i, u, n_img, n_cl, covis)
    g_rvec, trans, pts, posed = grown
    rot = Rotation.from_rotvec(g_rvec).as_matrix()
    ok = posed[obs_i] & ~np.isnan(pts[:, 0])[obs_c]
    f, rvec, tvec, p_ba, keep, res, inl = bundle_adjust(
        obs_c[ok],
        obs_i[ok],
        u[ok],
        rot,
        trans,
        pts,
        f_try,
        n_img,
        n_cl,
        opt_f=True,
    )
    full_keep = np.zeros(len(obs_c), bool)
    full_keep[np.nonzero(ok)[0]] = keep
    full_res = np.full(len(obs_c), np.inf)
    full_res[np.nonzero(ok)[0]] = res
    pts, keep, res = p_ba, full_keep, full_res
    rk = res[keep]
    n_pts = len(np.unique(obs_c[keep]))
    print(
        f"\nbootstrap result: f = {f:.1f} px, {n_pts} points, "
        f"{keep.sum()}/{len(obs_c)} observations kept"
    )
    print(
        f"reprojection (kept): rms {np.sqrt((rk**2).mean()):.2f} px, "
        f"median {np.median(rk):.2f} px; inlier<2px {100 * (res < 2).mean():.1f}% "
        f"of all obs"
    )

    compare_to_reference(data["names"], rvec, tvec, f)

    out = WS / "sfmr" / "bootstrap-pinhole.sfmr"
    save_sfmr(data, f, rvec, tvec, pts, keep, res, out)


if __name__ == "__main__":
    main()
