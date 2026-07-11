# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment E: coarse bootstrap from clusters + affine correspondences.

The question (see specs/core/grid-distortion-model.md): a single physical
camera captures a scene from many positions/angles; before any solve exists
we have track clusters and cluster-patch affine warps. Can they bootstrap
coarse cameras plus a coarse radial correction — no solver initialization
anywhere — well enough to be useful?

E1 (positions only): alternate a missing-data affine factorization (rank-3
per-image [M_i | t_i] / per-cluster X_c, ALS) with an L0 radial-profile fit
to the factorization residuals. Score against a reference solve of the same
subset: radial correction profile (modulo the linear/gauge term) and camera
row-space angles vs the solved rotations.

E2 (what the warps add): correct the observed member warps by the estimated
radial profile's Jacobian and measure whether the multi-cluster warp
factorization residual drops — i.e., do the affine correspondences sense
the recovered distortion? Plus a cluster-count sweep for E1.

Run: pixi run -e dev python scripts/exp_cluster_bootstrap.py <workspace>
(default workspace: dino_ledge_e_ws)
"""

import json
import sys
from pathlib import Path

import numpy as np
import scipy.interpolate as si

WS = Path(sys.argv[1] if len(sys.argv) > 1 else "dino_ledge_e_ws")
# Optional key=val args: rmax_frac=0.45 (cut observations beyond this
# fraction of the half-diagonal — fisheye corners diverge in the perspective
# chart), knots=6, init_f=<focal px> (seed the profile with the
# equidistant→perspective correction for that nominal focal — the production
# "base + delta from config" setting), min_span=6.
ARGS = dict(kv.split("=", 1) for kv in sys.argv[2:])
RMAX_FRAC = float(ARGS.get("rmax_frac", 1.0))
N_KNOTS = int(ARGS.get("knots", 6))
INIT_F = float(ARGS["init_f"]) if "init_f" in ARGS else None
MIN_SPAN = int(ARGS.get("min_span", 6))
RNG = np.random.default_rng(11)


# ── Small helpers ────────────────────────────────────────────────────────────


def design_1d(x, n_ctrl, lo, hi, degree=3):
    t = np.concatenate(
        [[lo] * degree, np.linspace(lo, hi, n_ctrl - degree + 1), [hi] * degree]
    )
    x = np.clip(x, lo, hi - 1e-9 * (hi - lo))
    return si.BSpline.design_matrix(x, t, degree).toarray()


def load_data():
    """Clusters + warps + keypoint positions + image dims from the workspace."""
    from sfmtool._sfmtool.io import read_matches, read_sift, read_sift_metadata

    patches = sorted((WS / "matches").glob("*-clusters-patches.matches"))
    assert patches, f"no cluster-patches file under {WS}/matches"
    data = read_matches(patches[0])
    names = list(data["image_names"])
    feature_counts = data["feature_counts"]
    prefix = data["metadata"]["workspace"]["contents"]["feature_prefix_dir"]

    positions = []
    dims = []
    for i, name in enumerate(names):
        rel = Path(name)
        sift_path = WS / rel.parent / prefix / f"{rel.name}.sift"
        meta = read_sift_metadata(sift_path)["metadata"]
        dims.append((meta["image_width"], meta["image_height"]))
        s = read_sift(sift_path)
        positions.append(
            np.ascontiguousarray(s["positions_xy"][: int(feature_counts[i])])
        )
    return data, names, positions, dims


def select_clusters(data, min_span):
    """Cluster member table for clusters with >= min_span usable members.

    Returns per-cluster arrays of (image_idx, feature_idx, is_reference,
    warp 2x2, zncc) using kept + reference members only.
    """
    starts = np.asarray(data["cluster_starts"])
    mi = np.asarray(data["member_images"])
    mf = np.asarray(data["member_features"])
    st = np.asarray(data["member_status"])
    aff = np.asarray(data["member_affines"])  # (M, 2, 3) absolute x_m = A x_ref + t
    zncc = np.asarray(data["member_zncc"])
    out = []
    for c in range(len(starts) - 1):
        lo, hi = starts[c], starts[c + 1]
        sel = np.nonzero((st[lo:hi] == 0) | (st[lo:hi] == 1))[0] + lo
        if len(sel) < min_span:
            continue
        out.append(
            {
                "img": mi[sel],
                "feat": mf[sel],
                "is_ref": st[sel] == 0,
                "A": aff[sel, :, :2],
                "zncc": zncc[sel],
            }
        )
    return out


# ── E1: alternating affine factorization + radial profile ───────────────────


def build_observations(clusters, positions, dims):
    """Per-observation (cluster, image, centered px position, radius norm)."""
    half = np.array([[w / 2, h / 2] for (w, h) in dims])
    rmax = np.linalg.norm(half[0])
    obs_c, obs_i, obs_u = [], [], []
    for c, cl in enumerate(clusters):
        for img, feat in zip(cl["img"], cl["feat"]):
            u = positions[img][feat] - half[img]
            if np.linalg.norm(u) > RMAX_FRAC * rmax:
                continue
            obs_c.append(c)
            obs_i.append(img)
            obs_u.append(u)
    return (
        np.asarray(obs_c),
        np.asarray(obs_i),
        np.asarray(obs_u, dtype=np.float64),
        rmax,
    )


def radial_correct(u, rmax, knots):
    """Apply the radial correction d(u) = g(r)·rhat to centered px positions."""
    r = np.linalg.norm(u, axis=1)
    rhat = np.divide(u, np.maximum(r, 1e-9)[:, None])
    des = design_1d(r / rmax, len(knots), 0.0, 1.0)
    return u + (des @ knots)[:, None] * rhat


def fit_radial(u, target_delta, rmax, n_knots):
    """Fit g(r) so that g(r)·rhat ≈ target_delta, then remove the r-linear
    component (the scale gauge the factorization absorbs)."""
    r = np.linalg.norm(u, axis=1)
    rhat = np.divide(u, np.maximum(r, 1e-9)[:, None])
    g_target = np.einsum("ij,ij->i", target_delta, rhat)
    des = design_1d(r / rmax, n_knots, 0.0, 1.0)
    knots = np.linalg.lstsq(des, g_target, rcond=None)[0]
    # Project out the linear-in-r component over the observed radii.
    gvals = des @ knots
    alpha = (gvals @ r) / (r @ r)
    knots -= np.linalg.lstsq(des, alpha * r, rcond=None)[0]
    return knots


def als_factorize(obs_c, obs_i, u_corr, n_img, n_cl, rounds=25, trim=0.0):
    """Missing-data affine factorization: u ≈ M_i X_c + t_i (M 2x3, X 3).

    ALS with per-image camera solves and per-cluster point solves. Init:
    X from the centered measurement SVD with mean-filling. Optionally trims
    the worst `trim` fraction of observations (robustness) in later rounds.
    Returns (M, t, X, per-obs residuals, inlier mask).
    """
    # Init X via mean-filled SVD.
    grid = np.zeros((2 * n_img, n_cl))
    cnt = np.zeros((n_img, n_cl))
    for c, i, u in zip(obs_c, obs_i, u_corr):
        grid[2 * i : 2 * i + 2, c] = u
        cnt[i, c] = 1
    row_mean = grid.sum(axis=1, keepdims=True) / np.maximum(
        np.repeat(cnt.sum(axis=1), 2)[:, None], 1
    )
    filled = np.where(np.repeat(cnt, 2, axis=0) > 0, grid, row_mean)
    filled -= row_mean
    _, _, vt = np.linalg.svd(filled, full_matrices=False)
    x_pts = vt[:3].T  # (n_cl, 3) initial points (arbitrary affine frame)

    keep = np.ones(len(obs_c), bool)
    m_cam = np.zeros((n_img, 2, 3))
    t_cam = np.zeros((n_img, 2))
    for it in range(rounds):
        # Cameras from points.
        for i in range(n_img):
            s = keep & (obs_i == i)
            if s.sum() < 4:
                continue
            xh = np.concatenate([x_pts[obs_c[s]], np.ones((s.sum(), 1))], axis=1)
            sol = np.linalg.lstsq(xh, u_corr[s], rcond=None)[0]  # (4, 2)
            m_cam[i] = sol[:3].T
            t_cam[i] = sol[3]
        # Points from cameras.
        for c in range(n_cl):
            s = keep & (obs_c == c)
            if s.sum() < 3:
                continue
            a = m_cam[obs_i[s]].reshape(-1, 3)
            b = (u_corr[s] - t_cam[obs_i[s]]).reshape(-1)
            x_pts[c] = np.linalg.lstsq(a, b, rcond=None)[0]
        res = u_corr - (
            np.einsum("nij,nj->ni", m_cam[obs_i], x_pts[obs_c]) + t_cam[obs_i]
        )
        if trim > 0 and it >= rounds // 2:
            rn = np.linalg.norm(res, axis=1)
            keep = rn < np.quantile(rn[keep], 1 - trim)
    return m_cam, t_cam, x_pts, res, keep


def bootstrap(obs_c, obs_i, obs_u, rmax, n_img, n_cl, n_knots=6, outer=4, knots0=None):
    """Alternate factorization and radial-profile fitting from scratch (or
    from a config-seeded initial profile)."""
    knots = np.zeros(n_knots) if knots0 is None else knots0.copy()
    for _ in range(outer):
        u_corr = radial_correct(obs_u, rmax, knots)
        m_cam, t_cam, x_pts, res, keep = als_factorize(
            obs_c, obs_i, u_corr, n_img, n_cl, trim=0.05
        )
        # New correction target: what delta would zero the residuals.
        target = radial_correct(obs_u, rmax, knots)[keep] + res[keep] - obs_u[keep]
        knots = fit_radial(obs_u[keep], target, rmax, n_knots)
    u_corr = radial_correct(obs_u, rmax, knots)
    m_cam, t_cam, x_pts, res, keep = als_factorize(
        obs_c, obs_i, u_corr, n_img, n_cl, trim=0.05
    )
    return knots, m_cam, res, keep


# ── Reference comparison ─────────────────────────────────────────────────────


def load_reference(names):
    """Reference camera + per-image rotations.

    Prefers the workspace's solved .sfmr; falls back to a
    `ref_camera.json` ({model, parameters}) with no rotations — used for
    the fisheye set, where the reference is the rig calibration rather
    than a solve."""
    from scipy.spatial.transform import Rotation

    from sfmtool._sfmtool import SfmrReconstruction

    sfmrs = sorted((WS / "sfmr").glob("*.sfmr"))
    if not sfmrs:
        ref = json.loads((WS / "ref_camera.json").read_text())
        return ref, [None] * len(names)
    rec = SfmrReconstruction.load(str(sfmrs[-1]))
    cam = rec.cameras[0]
    quats = np.asarray(rec.quaternions_wxyz)
    name_to_rot = {}
    for i in range(rec.image_count):
        q = quats[i]
        rot = Rotation.from_quat([q[1], q[2], q[3], q[0]]).as_matrix()
        name_to_rot[Path(rec.image_names[i]).name] = rot
    rots = [name_to_rot.get(Path(n).name) for n in names]
    return cam.to_dict(), rots


def reference_profile(cam, r_px):
    """Reference model's radial correction g_ref(r) in px (distorted→ideal
    perspective), from a {model, parameters} dict. SIMPLE_RADIAL and
    OPENCV_FISHEYE supported; returns (None, model) otherwise."""
    params = cam["parameters"]
    model = str(cam["model"]).upper()
    if "SIMPLE_RADIAL" in model:
        f = params["focal_length"]
        k = params["radial_distortion_k1"]
        # ideal x maps to distorted x(1 + k x^2) (normalized). Invert per r.
        rn_d = r_px / f
        rn_i = rn_d.copy()
        for _ in range(20):
            rn_i = rn_d / (1 + k * rn_i * rn_i)
        return (rn_i - rn_d) * f, model
    if "FISHEYE" in model:
        f = params["fx"]
        ks = np.array([params["k1"], params["k2"], params["k3"], params["k4"]])
        thetad = r_px / f
        theta = thetad.copy()
        for _ in range(30):  # invert thetad = theta(1 + k1 t^2 + ... )
            t2 = theta * theta
            poly = 1 + t2 * (ks[0] + t2 * (ks[1] + t2 * (ks[2] + t2 * ks[3])))
            dpoly = (
                3 * ks[0] * t2
                + 5 * ks[1] * t2**2
                + 7 * ks[2] * t2**3
                + 9 * ks[3] * t2**4
            )
            theta -= (theta * poly - thetad) / (poly + dpoly)
        return f * np.tan(theta) - r_px, model
    return None, model


def init_knots_equidistant(n_knots, rmax, f):
    """Seed profile: the equidistant→perspective correction for focal `f`
    (`g0(r) = f·tan(r/f) − r`), fitted on the observed radius range."""
    r = np.linspace(0.01, RMAX_FRAC, 200) * rmax
    theta = np.clip(r / f, 0, 1.45)
    g0 = f * np.tan(theta) - r
    des = design_1d(r / rmax, n_knots, 0.0, 1.0)
    return np.linalg.lstsq(des, g0, rcond=None)[0]


# ── E2: do the warps sense the recovered distortion? ─────────────────────────


def warp_residual(clusters, positions, dims, rmax, knots):
    """Multi-cluster warp factorization residual, with warps corrected by the
    current radial profile's Jacobian.

    Per cluster: members' absolute warps A (x_m = A x_ref) are corrected to
    A' = J(u_m)^{-1} A J(u_ref) using the profile's analytic Jacobian, then a
    rank-2 consistency fit per cluster (each corrected warp should factor
    through a common tangent frame) scores the set; we report the mean
    relative misfit. Lower = the warps agree better with that profile.
    """
    half = np.array([[w / 2, h / 2] for (w, h) in dims])

    def jac(u):
        # J = I + dg/du for d(u) = g(r) rhat, g spline in r/rmax.
        r = float(np.linalg.norm(u))
        if r < 1e-6:
            return np.eye(2)
        rhat = u / r
        des = design_1d(np.array([r / rmax]), len(knots), 0.0, 1.0)
        eps = 1e-3 * rmax
        des2 = design_1d(np.array([(r + eps) / rmax]), len(knots), 0.0, 1.0)
        g = (des @ knots)[0]
        dg = ((des2 @ knots)[0] - g) / eps
        outer_rr = np.outer(rhat, rhat)
        return np.eye(2) + dg * outer_rr + (g / r) * (np.eye(2) - outer_rr)

    misfits = []
    for cl in clusters:
        refs = np.nonzero(cl["is_ref"])[0]
        if len(refs) != 1 or len(cl["img"]) < 4:
            continue
        u_ref = (
            positions[cl["img"][refs[0]]][cl["feat"][refs[0]]]
            - half[cl["img"][refs[0]]]
        )
        j_ref = jac(u_ref)
        mats = []
        for k_m in range(len(cl["img"])):
            if cl["is_ref"][k_m]:
                continue
            u_m = positions[cl["img"][k_m]][cl["feat"][k_m]] - half[cl["img"][k_m]]
            a_corr = np.linalg.solve(jac(u_m), cl["A"][k_m] @ j_ref)
            mats.append(a_corr.ravel())
        if len(mats) < 3:
            continue
        m = np.asarray(mats)  # (n, 4)
        # Rank-2 misfit of the stacked warps (weak-perspective consistency).
        _, s, _ = np.linalg.svd(m - m.mean(0), full_matrices=False)
        denom = np.linalg.norm(m - m.mean(0))
        if denom > 1e-9:
            misfits.append(np.sqrt((s[2:] ** 2).sum()) / denom)
    return float(np.mean(misfits)), len(misfits)


# ── Main ─────────────────────────────────────────────────────────────────────


def main():
    data, names, positions, dims = load_data()
    n_img = len(names)
    clusters_all = select_clusters(data, min_span=MIN_SPAN)
    print(
        f"{WS}: {n_img} images; {len(clusters_all)} clusters with span >= "
        f"{MIN_SPAN}; rmax_frac={RMAX_FRAC} knots={N_KNOTS} init_f={INIT_F}"
    )

    cam, rots = load_reference(names)
    print(f"reference camera: {cam}")

    for n_cl_target in [150, 500, len(clusters_all)]:
        if n_cl_target > len(clusters_all):
            continue
        idx = RNG.choice(
            len(clusters_all), min(n_cl_target, len(clusters_all)), replace=False
        )
        clusters = [clusters_all[i] for i in sorted(idx)]
        obs_c, obs_i, obs_u, rmax = build_observations(clusters, positions, dims)
        knots0 = init_knots_equidistant(N_KNOTS, rmax, INIT_F) if INIT_F else None
        knots, m_cam, res, keep = bootstrap(
            obs_c,
            obs_i,
            obs_u,
            rmax,
            n_img,
            len(clusters),
            n_knots=N_KNOTS,
            knots0=knots0,
        )
        rn = np.linalg.norm(res[keep], axis=1)

        # Score the radial profile vs the reference model, modulo linear.
        r_eval = np.linspace(0.15, 0.95, 60) * RMAX_FRAC * rmax
        des = design_1d(r_eval / rmax, len(knots), 0.0, 1.0)
        g_est = des @ knots
        g_ref, model = reference_profile(cam, r_eval)
        if g_ref is not None:
            # Remove each profile's linear component before comparing.
            for g in (g_est, g_ref):
                g -= (g @ r_eval) / (r_eval @ r_eval) * r_eval
            prof_rms = np.sqrt(((g_est - g_ref) ** 2).mean())
            prof_scale = np.sqrt((g_ref**2).mean())
        else:
            prof_rms, prof_scale = float("nan"), float("nan")

        # Camera comparison vs reference rotations. The factorization's
        # cameras carry a global 3x3 affine gauge (M_i' = M_i G undoes any
        # G applied to the points), so first solve one G plus per-image
        # scales aligning M_i G ~ s_i R_i[:2] (weak-perspective form), then
        # measure the per-camera row-space angle. Both quaternion
        # conventions are tried; the better mean is reported.
        def gauge_aligned_angles(rows_of):
            valid = [
                i
                for i in range(n_img)
                if rots[i] is not None and np.linalg.norm(m_cam[i]) > 1e-9
            ]
            m_stack = np.concatenate([m_cam[i] for i in valid])
            s_i = {i: 1.0 for i in valid}
            g = np.eye(3)
            for _ in range(5):
                r_stack = np.concatenate([s_i[i] * rows_of(i) for i in valid])
                g = np.linalg.lstsq(m_stack, r_stack, rcond=None)[0]
                for i in valid:
                    mg = m_cam[i] @ g
                    denom = (mg * mg).sum()
                    if denom > 1e-12:
                        s_i[i] = (mg * rows_of(i)).sum() / denom
            out = []
            for i in valid:
                qm, _ = np.linalg.qr((m_cam[i] @ g).T)
                qr_, _ = np.linalg.qr(rows_of(i).T)
                s = np.linalg.svd(qm.T @ qr_, compute_uv=False)
                out.append(np.degrees(np.arccos(np.clip(s.min(), -1, 1))))
            return np.mean(out)

        if any(r is not None for r in rots):
            angle_mean = min(
                gauge_aligned_angles(lambda i: np.asarray(rots[i])[:2]),
                gauge_aligned_angles(lambda i: np.asarray(rots[i]).T[:2]),
            )
        else:
            angle_mean = float("nan")

        # E2: warp agreement with zero vs estimated profile.
        mis0, n0 = warp_residual(clusters, positions, dims, rmax, np.zeros_like(knots))
        mis1, _ = warp_residual(clusters, positions, dims, rmax, knots)

        print(
            f"clusters={len(clusters):5d} obs={len(obs_c):6d} | factorization "
            f"residual rms {np.sqrt((rn**2).mean()):6.2f} px | profile err vs "
            f"{model} (mod linear) {prof_rms:6.2f} px (ref profile scale "
            f"{prof_scale:.2f} px) | cam row-space angle mean "
            f"{angle_mean:5.2f} deg | warp rank-2 misfit: zero-profile "
            f"{mis0:.4f} vs est-profile {mis1:.4f} ({n0} clusters)"
        )


if __name__ == "__main__":
    main()
