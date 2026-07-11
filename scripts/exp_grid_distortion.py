# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Grid-distortion model experiments (see specs/core/grid-distortion-model.md).

Experiment A — expressiveness: fit the hierarchical grid (L0 radial profile,
L1 quadrant-symmetric spline, L2 full spline; knot doubling within levels)
against known parametric distortion fields and report DOF vs residual.

Experiment B — absorption/degeneracy: a synthetic bundle adjustment with
poses + points + focal + a full grid, under full vs center-only coverage and
with the linear-gauge constraint on/off; reports structure / focal / field
recovery error per condition.

Run: pixi run python scripts/exp_grid_distortion.py
"""

import json
from pathlib import Path

import numpy as np
import scipy.interpolate as si
import scipy.optimize
import scipy.sparse
from scipy.spatial.transform import Rotation

RNG = np.random.default_rng(7)


# ── B-spline helpers ─────────────────────────────────────────────────────────


def open_knots(n_ctrl: int, lo: float, hi: float, degree: int = 3) -> np.ndarray:
    """Open-uniform knot vector for `n_ctrl` control points on [lo, hi]."""
    n_seg = n_ctrl - degree
    assert n_seg >= 1, "need at least degree+1 control points"
    inner = np.linspace(lo, hi, n_seg + 1)
    return np.concatenate([[lo] * degree, inner, [hi] * degree])


def design_1d(x: np.ndarray, n_ctrl: int, lo: float, hi: float, degree: int = 3):
    """Dense (len(x), n_ctrl) cubic B-spline design matrix on [lo, hi]."""
    t = open_knots(n_ctrl, lo, hi, degree)
    x = np.clip(x, lo, hi - 1e-9 * (hi - lo))
    return si.BSpline.design_matrix(x, t, degree).toarray()


def design_2d(xy: np.ndarray, n_ctrl: int, lo: float, hi: float) -> np.ndarray:
    """Dense (N, n_ctrl^2) tensor-product cubic design matrix on [lo,hi]^2.

    Column order matches `coef.reshape(n_ctrl, n_ctrl)` with x-major rows
    (kron(row_x, row_y))."""
    bx = design_1d(xy[:, 0], n_ctrl, lo, hi)
    by = design_1d(xy[:, 1], n_ctrl, lo, hi)
    return np.einsum("ni,nj->nij", bx, by).reshape(xy.shape[0], -1)


def lstsq_field(design: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Least-squares fit of a scalar field; returns the fitted values."""
    coef = np.linalg.lstsq(design, target, rcond=None)[0]
    return design @ coef


# ── Experiment A: expressiveness vs known models ─────────────────────────────


def fisheye_theta_from_thetad(thetad: np.ndarray, ks: np.ndarray) -> np.ndarray:
    """Invert OPENCV_FISHEYE's thetad = theta(1 + k1 t^2 + ... + k4 t^8)."""
    theta = thetad.copy()
    for _ in range(30):
        t2 = theta * theta
        poly = 1 + t2 * (ks[0] + t2 * (ks[1] + t2 * (ks[2] + t2 * ks[3])))
        dpoly = (
            3 * ks[0] * t2 + 5 * ks[1] * t2**2 + 7 * ks[2] * t2**3 + 9 * ks[3] * t2**4
        )
        f = theta * poly - thetad
        theta = theta - f / (poly + dpoly)
    return theta


def field_fisheye_kerry(n: int = 91):
    """Kerry Park OPENCV_FISHEYE delta field in the azimuthal angle chart.

    Samples observed pixels inside the image circle; the field is
    `(theta - thetad) * rhat`, i.e. the correction the grid must carry on
    top of the equidistant base, scaled by f to pixels.
    """
    cfg = json.loads(
        (
            Path(__file__).parent.parent / "test-data/images/kerry_park/rig_config.json"
        ).read_text()
    )
    cam = cfg[0]["cameras"][0]
    fx, fy, cx, cy, *ks = cam["camera_params"]
    ks = np.asarray(ks)
    f = 0.5 * (fx + fy)

    u = np.linspace(4, 476, n)
    uu, vv = np.meshgrid(u, u)
    du, dv = (uu - cx) / f, (vv - cy) / f
    rd = np.hypot(du, dv)
    keep = rd < (236.0 / f)
    du, dv, rd = du[keep], dv[keep], rd[keep]
    thetad = rd
    theta = fisheye_theta_from_thetad(thetad, ks)
    scale = np.where(rd > 1e-12, (theta - thetad) / np.maximum(rd, 1e-12), 0.0)
    dx = scale * du * f  # px-equivalent (angle error * f)
    dy = scale * dv * f
    # Chart coordinates in [-1, 1] over the image circle.
    lim = 236.0 / f
    xy = np.stack([du / lim, dv / lim], axis=1)
    return "kerry_fisheye", xy, np.stack([dx, dy], axis=1)


def field_pixel_model(name: str, dist_fn, f: float, w: float, h: float, n: int = 91):
    """Displacement field (px) of a narrow-FOV model over the image plane."""
    xs = np.linspace(-w / 2, w / 2, n) / f
    ys = np.linspace(-h / 2, h / 2, n) / f
    xx, yy = np.meshgrid(xs, ys)
    x, y = xx.ravel(), yy.ravel()
    dxn, dyn = dist_fn(x, y)
    lim = max(w, h) / 2 / f
    xy = np.stack([x / lim, y / lim], axis=1)
    return name, xy, np.stack([(dxn - x) * f, (dyn - y) * f], axis=1)


def simple_radial(k1):
    def fn(x, y):
        r2 = x * x + y * y
        s = 1 + k1 * r2
        return x * s, y * s

    return fn


def brown_tangential(k1, k2, p1, p2):
    def fn(x, y):
        r2 = x * x + y * y
        s = 1 + k1 * r2 + k2 * r2 * r2
        xt = 2 * p1 * x * y + p2 * (r2 + 2 * x * x)
        yt = p1 * (r2 + 2 * y * y) + 2 * p2 * x * y
        return x * s + xt, y * s + yt

    return fn


def symmetrize(xy: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Project the field onto its 4-fold covariant-symmetric part.

    Averages the field over the four axis reflections with vector parity
    (`d(−x, y) = (−dx, dy)` etc.), interpolating the reflected samples from
    the original grid (samples come from a symmetric meshgrid, so reflected
    sample positions are exact grid positions; use nearest lookup by
    rounding through a dict for robustness).
    """
    key = {(round(p[0], 9), round(p[1], 9)): i for i, p in enumerate(xy)}
    out = np.zeros_like(d)
    n_found = 0
    for i, (x, y) in enumerate(xy):
        acc = d[i].copy()
        cnt = 1
        for sx, sy in [(-1, 1), (1, -1), (-1, -1)]:
            j = key.get((round(sx * x, 9), round(sy * y, 9)))
            if j is not None:
                acc[0] += sx * d[j, 0]
                acc[1] += sy * d[j, 1]
                cnt += 1
        out[i] = acc / cnt
        n_found += cnt == 4
    return out


def fit_level(xy, d, level: str, n_ctrl: int):
    """Fit one hierarchy level; returns (dof, fitted field)."""
    if level == "L0":  # radial profile g(r), d = g(r) * rhat
        r = np.linalg.norm(xy, axis=1)
        rhat = np.divide(xy, np.maximum(r, 1e-12)[:, None])
        g = np.where(r > 1e-12, np.einsum("ij,ij->i", d, rhat), 0.0)
        des = design_1d(r, n_ctrl, 0.0, float(r.max()))
        gfit = lstsq_field(des, g)
        return n_ctrl, gfit[:, None] * rhat
    if level == "L1":  # quadrant-symmetric: fit the symmetrized field
        dsym = symmetrize(xy, d)
        des = design_2d(xy, n_ctrl, -1.0, 1.0)
        fit = np.stack(
            [lstsq_field(des, dsym[:, 0]), lstsq_field(des, dsym[:, 1])], axis=1
        )
        # Quadrant tying: n_ctrl^2 coefficients per component, ~1/4 free.
        dof = 2 * int(np.ceil(n_ctrl / 2)) ** 2
        return dof, fit
    if level == "L2":
        des = design_2d(xy, n_ctrl, -1.0, 1.0)
        fit = np.stack([lstsq_field(des, d[:, 0]), lstsq_field(des, d[:, 1])], axis=1)
        return 2 * n_ctrl * n_ctrl, fit
    raise ValueError(level)


def experiment_a():
    print("=" * 72)
    print("Experiment A: hierarchy expressiveness vs known parametric models")
    print("=" * 72)
    models = [
        field_fisheye_kerry(),
        field_pixel_model(
            "simple_radial(k=-0.12)", simple_radial(-0.12), 600, 800, 600
        ),
        field_pixel_model(
            "brown+tangential(k1=-.28,k2=.07,p1=.0012,p2=-.0007)",
            brown_tangential(-0.28, 0.07, 0.0012, -0.0007),
            600,
            800,
            600,
        ),
    ]
    for name, xy, d in models:
        mag = np.linalg.norm(d, axis=1)
        anti = d - symmetrize(xy, d)
        print(f"\n--- {name}")
        print(
            f"    field magnitude: rms {np.sqrt((mag**2).mean()):8.3f} px, "
            f"max {mag.max():8.3f} px"
        )
        anti_mag = np.linalg.norm(anti, axis=1)
        print(
            f"    antisymmetric part (L1 floor): rms "
            f"{np.sqrt((anti_mag**2).mean()):.4f} px, max {anti_mag.max():.4f} px"
        )
        print(
            f"    {'level':6s} {'ctrl':>5s} {'DOF':>5s} {'rms_px':>9s} {'max_px':>9s}"
        )
        for level, ctrls in [
            ("L0", [4, 6, 8, 12]),
            ("L1", [4, 6, 8, 12]),
            ("L2", [4, 6, 8, 12]),
        ]:
            for c in ctrls:
                dof, fit = fit_level(xy, d, level, c)
                err = np.linalg.norm(fit - d, axis=1)
                print(
                    f"    {level:6s} {c:5d} {dof:5d} "
                    f"{np.sqrt((err**2).mean()):9.4f} {err.max():9.4f}"
                )


# ── Experiment B: absorption / degeneracy in a mock BA ──────────────────────

W, H, F_GT = 800.0, 600.0, 600.0
N_CTRL_B = 6  # full-grid control points per axis (2*36 = 72 grid DOF)


def grid_displace(xy_norm: np.ndarray, coef: np.ndarray) -> np.ndarray:
    """Fast tensor-spline evaluation: (Bx @ C) row-dotted with By."""
    bx = design_1d(xy_norm[:, 0], N_CTRL_B, -1.0, 1.0)
    by = design_1d(xy_norm[:, 1], N_CTRL_B, -1.0, 1.0)
    cx = coef[: N_CTRL_B**2].reshape(N_CTRL_B, N_CTRL_B)
    cy = coef[N_CTRL_B**2 :].reshape(N_CTRL_B, N_CTRL_B)
    return np.stack(
        [((bx @ cx) * by).sum(axis=1), ((bx @ cy) * by).sum(axis=1)], axis=1
    )


def make_scene(mode: str, coverage: str):
    """Synthetic scene.

    `frontal`: points in a slab, cameras translating/nodding in front of it —
    the weak-geometry regime where radial distortion and depth trade off (the
    "dome" degeneracy). `orbit`: cameras on a ring looking inward at a point
    ball — a real object-sweep capture, geometry that should pin distortion.
    """
    n_pts, n_cams = 400, 20
    poses = []
    if mode == "frontal":
        pts = np.stack(
            [
                RNG.uniform(-4, 4, n_pts),
                RNG.uniform(-3, 3, n_pts),
                RNG.uniform(8, 14, n_pts),
            ],
            axis=1,
        )
        for i in range(n_cams):
            ang = 0.25 * np.sin(2 * np.pi * i / n_cams)
            rot = Rotation.from_euler(
                "yxz", [ang, 0.08 * np.cos(2 * np.pi * i / n_cams), 0]
            )
            t = np.array(
                [
                    2.5 * np.sin(2 * np.pi * i / n_cams),
                    0.8 * np.cos(4 * np.pi * i / n_cams),
                    0,
                ]
            )
            poses.append((rot.as_rotvec(), t))
    elif mode == "orbit":
        pts = RNG.uniform(-3, 3, (n_pts, 3))
        for i in range(n_cams):
            a = 2 * np.pi * i / n_cams
            center = np.array([10 * np.cos(a), 2.0 * np.sin(2 * a), 10 * np.sin(a)])
            fwd = -center / np.linalg.norm(center)
            right = np.cross([0.0, 1.0, 0.0], fwd)
            right /= np.linalg.norm(right)
            up = np.cross(fwd, right)
            r_mat = np.stack([right, up, fwd])  # world -> camera rows
            rvec = Rotation.from_matrix(r_mat).as_rotvec()
            tvec = -r_mat @ center
            poses.append((rvec, tvec))
    else:
        raise ValueError(mode)
    lim = 0.45 if coverage == "center" else 1.0
    return pts, poses, lim


def project(rvec, tvec, pts, f, grid_coef, dist_fn=None):
    """Pinhole + (grid or parametric) distortion, in pixels from center."""
    q = Rotation.from_rotvec(rvec).apply(pts) + tvec
    x, y = q[:, 0] / q[:, 2], q[:, 1] / q[:, 2]
    if dist_fn is not None:
        x, y = dist_fn(x, y)
    u = np.stack([f * x, f * y], axis=1)
    if grid_coef is not None:
        xy_norm = u / np.array([W / 2, H / 2])
        u = u + grid_displace(np.clip(xy_norm, -1, 1), grid_coef)
    return u


def map_error(got: np.ndarray, gt: np.ndarray):
    """Camera-map residual modulo one uniform scale.

    Monocular BA has a global f↔scene-scale gauge (scaling f while the
    structure rescales is reprojection-invariant and removed from the
    structure metric by the similarity alignment), so map recovery must be
    judged modulo the best uniform scale between the learned and reference
    maps. Returns (per-sample error, fitted scale)."""
    a = float((got * gt).sum() / (got * got).sum())
    return np.linalg.norm(a * got - gt, axis=1), a


def experiment_b():
    print("\n" + "=" * 72)
    print("Experiment B: BA absorption/degeneracy (poses + points + f + full grid)")
    print("=" * 72)
    gt_dist = simple_radial(-0.10)
    noise = 0.3

    conditions = [
        ("frontal", "full", True),
        ("frontal", "full", False),
        ("frontal", "center", True),
        ("orbit", "full", True),
        ("orbit", "full", False),
        ("orbit", "center", True),
    ]
    for mode, coverage, gauge in conditions:
        if True:
            pts_gt, poses_gt, lim = make_scene(mode, coverage)
            # Observations: GT parametric distortion, no grid. Also record each
            # observation's ideal (pre-distortion) normalized coordinates —
            # the domain the composite-map recovery metric evaluates on.
            obs, vis, ideal = [], [], []
            for ci, (rv, tv) in enumerate(poses_gt):
                q = Rotation.from_rotvec(rv).apply(pts_gt) + tv
                xn, yn = q[:, 0] / q[:, 2], q[:, 1] / q[:, 2]
                u = project(rv, tv, pts_gt, F_GT, None, gt_dist)
                inside = (np.abs(u[:, 0]) < lim * W / 2) & (
                    np.abs(u[:, 1]) < lim * H / 2
                )
                for pi in np.nonzero(inside)[0]:
                    obs.append(u[pi] + RNG.normal(0, noise, 2))
                    vis.append((ci, pi))
                    ideal.append((xn[pi], yn[pi]))
            obs = np.asarray(obs)
            vis = np.asarray(vis)
            ideal = np.asarray(ideal)

            n_cams, n_pts = len(poses_gt), len(pts_gt)
            n_grid = 2 * N_CTRL_B**2

            def pack(poses, pts, f, grid):
                return np.concatenate(
                    [
                        np.concatenate([np.concatenate(p) for p in poses]),
                        pts.ravel(),
                        [f],
                        grid,
                    ]
                )

            def unpack(x):
                p = x[: 6 * n_cams].reshape(n_cams, 6)
                pts = x[6 * n_cams : 6 * n_cams + 3 * n_pts].reshape(n_pts, 3)
                f = x[6 * n_cams + 3 * n_pts]
                grid = x[6 * n_cams + 3 * n_pts + 1 :]
                return p, pts, f, grid

            # Gauge residuals: penalize the grid's best-fit constant + linear
            # component (the modes indistinguishable from f/PP/skew).
            xy_lin = np.stack(
                np.meshgrid(np.linspace(-1, 1, 9), np.linspace(-1, 1, 9)), axis=-1
            ).reshape(-1, 2)
            des_lin = design_2d(xy_lin, N_CTRL_B, -1, 1)
            basis = np.stack([np.ones(len(xy_lin)), xy_lin[:, 0], xy_lin[:, 1]], axis=1)
            # (3, n_ctrl^2): each grid coefficient's constant/x/y linear-fit
            # component — the K-indistinguishable modes the gauge pins to 0.
            proj_lin = np.linalg.lstsq(basis, des_lin, rcond=None)[0]

            def residuals(x):
                p, pts, f, grid = unpack(x)
                res = np.empty((len(obs), 2))
                for ci in range(n_cams):
                    sel = vis[:, 0] == ci
                    u = project(p[ci, :3], p[ci, 3:], pts[vis[sel, 1]], f, grid)
                    res[sel] = u - obs[sel]
                out = [res.ravel()]
                if gauge:
                    lin_x = proj_lin @ grid[: N_CTRL_B**2]
                    lin_y = proj_lin @ grid[N_CTRL_B**2 :]
                    out.append(50.0 * np.concatenate([lin_x, lin_y]))
                return np.concatenate(out)

            # Start from noisy GT poses/points, f off by 2%, zero grid.
            x0 = pack(
                [
                    (rv + RNG.normal(0, 0.002, 3), tv + RNG.normal(0, 0.02, 3))
                    for rv, tv in poses_gt
                ],
                pts_gt + RNG.normal(0, 0.03, pts_gt.shape),
                F_GT * 1.02,
                np.zeros(n_grid),
            )

            # Finite-difference sparsity: each residual touches its camera,
            # its point, f, and the grid.
            n_par = len(x0)
            rows, cols = [], []
            for oi, (ci, pi) in enumerate(vis):
                for r in (2 * oi, 2 * oi + 1):
                    cs = (
                        list(range(6 * ci, 6 * ci + 6))
                        + list(range(6 * n_cams + 3 * pi, 6 * n_cams + 3 * pi + 3))
                        + [6 * n_cams + 3 * n_pts]
                        + list(range(6 * n_cams + 3 * n_pts + 1, n_par))
                    )
                    rows += [r] * len(cs)
                    cols += cs
            n_res = 2 * len(obs) + (6 if gauge else 0)
            if gauge:
                for r in range(2 * len(obs), n_res):
                    rows += [r] * n_grid
                    cols += list(range(6 * n_cams + 3 * n_pts + 1, n_par))
            spar = scipy.sparse.coo_matrix(
                (np.ones(len(rows)), (rows, cols)), shape=(n_res, n_par)
            )

            sol = scipy.optimize.least_squares(
                residuals,
                x0,
                jac_sparsity=spar,
                tr_solver="lsmr",
                xtol=1e-10,
                max_nfev=60,
            )
            p, pts, f, grid = unpack(sol.x)

            # Structure error after similarity alignment to GT.
            mu_a, mu_b = pts.mean(0), pts_gt.mean(0)
            a, b = pts - mu_a, pts_gt - mu_b
            s = np.sqrt((b**2).sum() / (a**2).sum())
            u_, _, vt = np.linalg.svd(a.T @ b)
            rot = (u_ @ vt).T
            struct_rmse = np.sqrt(((s * a @ rot.T - b) ** 2).sum(axis=1).mean())

            # Camera-map recovery: compare the learned composite projection
            # `f·x + grid(f·x)` against the GT `F_GT·dist(x)` at the ideal
            # coordinates of the actual observations. Gauge-invariant (the
            # f↔grid-linear split cancels in the composite) and evaluated
            # only where data existed.
            sub = ideal[:: max(1, len(ideal) // 4000)]
            dxn, dyn = gt_dist(sub[:, 0], sub[:, 1])
            gt_map = np.stack([dxn, dyn], axis=1) * F_GT
            u_lin = sub * f
            got_map = u_lin + grid_displace(
                np.clip(u_lin / np.array([W / 2, H / 2]), -1, 1), grid
            )
            merr, _ = map_error(got_map, gt_map)

            rms_reproj = np.sqrt((sol.fun[: 2 * len(obs)] ** 2).mean())
            print(
                f"{mode:7s} coverage={coverage:6s} gauge={str(gauge):5s} | "
                f"reproj {rms_reproj:6.3f} px | "
                f"f {f:7.2f} (gt {F_GT:.0f}; carries the linear gauge) | "
                f"struct rmse {struct_rmse:7.4f} | camera-map err rms "
                f"{np.sqrt((merr**2).mean()):7.3f} max {merr.max():7.3f} px"
            )


def experiment_c():
    """Center-out radial curriculum on the absorbing (frontal) scene.

    Stage 1 solves poses + points + f from central observations only, grid
    frozen (all models — including 'none' — agree near the center, so this
    anchors structure with minimal distortion coupling). Later stages admit
    observations at growing radii and unlock the grid. Compared against the
    joint one-shot solve on the identical scene/observations.
    """
    print("\n" + "=" * 72)
    print("Experiment C: center-out curriculum vs joint solve (frontal, gauge on)")
    print("=" * 72)
    gt_dist = simple_radial(-0.10)
    noise = 0.3
    pts_gt, poses_gt, lim = make_scene("frontal", "full")
    n_cams, n_pts = len(poses_gt), len(pts_gt)
    n_grid = 2 * N_CTRL_B**2
    base = 6 * n_cams + 3 * n_pts + 1

    obs, vis, ideal = [], [], []
    for ci, (rv, tv) in enumerate(poses_gt):
        q = Rotation.from_rotvec(rv).apply(pts_gt) + tv
        xn, yn = q[:, 0] / q[:, 2], q[:, 1] / q[:, 2]
        u = project(rv, tv, pts_gt, F_GT, None, gt_dist)
        inside = (np.abs(u[:, 0]) < lim * W / 2) & (np.abs(u[:, 1]) < lim * H / 2)
        for pi in np.nonzero(inside)[0]:
            obs.append(u[pi] + RNG.normal(0, noise, 2))
            vis.append((ci, pi))
            ideal.append((xn[pi], yn[pi]))
    obs, vis, ideal = np.asarray(obs), np.asarray(vis), np.asarray(ideal)
    # Box-normalized radius of each observation (1.0 = image border).
    obs_r = np.maximum(np.abs(obs[:, 0]) / (W / 2), np.abs(obs[:, 1]) / (H / 2))

    def unpack(x):
        p = x[: 6 * n_cams].reshape(n_cams, 6)
        pts = x[6 * n_cams : 6 * n_cams + 3 * n_pts].reshape(n_pts, 3)
        return p, pts, x[base - 1], x[base:]

    xy_lin = np.stack(
        np.meshgrid(np.linspace(-1, 1, 9), np.linspace(-1, 1, 9)), axis=-1
    ).reshape(-1, 2)
    des_lin = design_2d(xy_lin, N_CTRL_B, -1, 1)
    lin_basis = np.stack([np.ones(len(xy_lin)), xy_lin[:, 0], xy_lin[:, 1]], axis=1)
    proj_lin = np.linalg.lstsq(lin_basis, des_lin, rcond=None)[0]

    def solve(x_start, sel, opt_grid, max_nfev):
        obs_s, vis_s = obs[sel], vis[sel]
        theta0 = x_start if opt_grid else x_start[:base]

        def to_full(theta):
            return theta if opt_grid else np.concatenate([theta, x_start[base:]])

        def residuals(theta):
            x = to_full(theta)
            p, pts, f, grid = unpack(x)
            res = np.empty((len(obs_s), 2))
            for ci in range(n_cams):
                m = vis_s[:, 0] == ci
                if m.any():
                    res[m] = (
                        project(p[ci, :3], p[ci, 3:], pts[vis_s[m, 1]], f, grid)
                        - obs_s[m]
                    )
            out = [res.ravel()]
            if opt_grid:
                out.append(
                    50.0
                    * np.concatenate(
                        [proj_lin @ grid[: N_CTRL_B**2], proj_lin @ grid[N_CTRL_B**2 :]]
                    )
                )
            return np.concatenate(out)

        n_par = len(theta0)
        rows, cols = [], []
        for oi, (ci, pi) in enumerate(vis_s):
            cs = (
                list(range(6 * ci, 6 * ci + 6))
                + list(range(6 * n_cams + 3 * pi, 6 * n_cams + 3 * pi + 3))
                + [base - 1]
                + (list(range(base, n_par)) if opt_grid else [])
            )
            for r in (2 * oi, 2 * oi + 1):
                rows += [r] * len(cs)
                cols += cs
        n_res = 2 * len(obs_s) + (6 if opt_grid else 0)
        if opt_grid:
            for r in range(2 * len(obs_s), n_res):
                rows += [r] * n_grid
                cols += list(range(base, n_par))
        spar = scipy.sparse.coo_matrix(
            (np.ones(len(rows)), (rows, cols)), shape=(n_res, n_par)
        )
        sol = scipy.optimize.least_squares(
            residuals,
            theta0,
            jac_sparsity=spar,
            tr_solver="lsmr",
            xtol=1e-10,
            max_nfev=max_nfev,
        )
        return to_full(sol.x)

    def report(label, x):
        p, pts, f, grid = unpack(x)
        mu_a, mu_b = pts.mean(0), pts_gt.mean(0)
        a, b = pts - mu_a, pts_gt - mu_b
        s = np.sqrt((b**2).sum() / (a**2).sum())
        u_, _, vt = np.linalg.svd(a.T @ b)
        rot = (u_ @ vt).T
        struct = np.sqrt(((s * a @ rot.T - b) ** 2).sum(axis=1).mean())
        sub = ideal[:: max(1, len(ideal) // 4000)]
        dxn, dyn = gt_dist(sub[:, 0], sub[:, 1])
        gt_map = np.stack([dxn, dyn], axis=1) * F_GT
        u_lin = sub * f
        got = u_lin + grid_displace(
            np.clip(u_lin / np.array([W / 2, H / 2]), -1, 1), grid
        )
        merr, _ = map_error(got, gt_map)
        print(
            f"{label:26s} | f {f:7.2f} | struct rmse {struct:7.4f} | "
            f"camera-map err rms {np.sqrt((merr**2).mean()):7.3f} "
            f"max {merr.max():7.3f} px"
        )

    x0 = np.concatenate(
        [
            np.concatenate(
                [
                    np.concatenate(
                        [rv + RNG.normal(0, 0.002, 3), tv + RNG.normal(0, 0.02, 3)]
                    )
                    for rv, tv in poses_gt
                ]
            ),
            (pts_gt + RNG.normal(0, 0.03, pts_gt.shape)).ravel(),
            [F_GT * 1.02],
            np.zeros(n_grid),
        ]
    )

    report("joint (one shot)", solve(x0, np.ones(len(obs), bool), True, 60))

    x = solve(x0, obs_r < 0.35, False, 40)
    report("  stage 1: r<0.35, no grid", x)
    x = solve(x, obs_r < 0.65, True, 40)
    report("  stage 2: r<0.65, +grid", x)
    x = solve(x, np.ones(len(obs), bool), True, 40)
    report("center-out (stage 3: all)", x)


SOUTH_BUILDING = Path("C:/DataSets/ColmapDataSets/south-building")


def experiment_d():
    """Real data: recover the COLMAP South Building camera with the grid.

    Loads the dataset's solved sparse model (one shared SIMPLE_RADIAL camera,
    128 images walking around a building — orbit-class geometry), takes a
    subset of images and well-tracked points with their *real* observed
    keypoints, and re-solves poses + points + f + full grid from a pinhole
    start (f off by 2%, grid zero). The COLMAP-solved model is the reference
    map the composite `f·x + grid` should recover.
    """
    print("\n" + "=" * 72)
    print("Experiment D: real data (COLMAP South Building), grid vs solved model")
    print("=" * 72)
    if not (SOUTH_BUILDING / "sparse/cameras.txt").exists():
        print(f"skipped: {SOUTH_BUILDING} not found")
        return
    import pycolmap

    global W, H
    w_save, h_save = W, H

    rec = pycolmap.Reconstruction(str(SOUTH_BUILDING / "sparse"))
    cam = next(iter(rec.cameras.values()))
    f_ref, cx, cy, k_ref = (float(p) for p in cam.params)
    W, H = float(cam.width), float(cam.height)
    ref_dist = simple_radial(k_ref)

    def ref_map(x):
        dxn, dyn = ref_dist(x[:, 0], x[:, 1])
        return np.stack([dxn, dyn], axis=1) * f_ref

    # Subset: 24 evenly spaced images; points tracked >= 4x among them.
    images = sorted(rec.images.values(), key=lambda im: im.name)
    sel_imgs = images[:: max(1, len(images) // 24)][:24]
    counts: dict[int, int] = {}
    for im in sel_imgs:
        for p2 in im.points2D:
            if p2.has_point3D():
                counts[p2.point3D_id] = counts.get(p2.point3D_id, 0) + 1
    kept_ids = sorted(pid for pid, c in counts.items() if c >= 4)[:800]
    pid_to_local = {pid: i for i, pid in enumerate(kept_ids)}
    pts0 = np.array([rec.points3D[pid].xyz for pid in kept_ids])

    poses, obs, vis = [], [], []
    for ci, im in enumerate(sel_imgs):
        m = np.asarray(im.cam_from_world().matrix())
        poses.append((Rotation.from_matrix(m[:, :3]).as_rotvec(), m[:, 3]))
        for p2 in im.points2D:
            if p2.has_point3D() and p2.point3D_id in pid_to_local:
                obs.append(np.asarray(p2.xy) - np.array([cx, cy]))
                vis.append((ci, pid_to_local[p2.point3D_id]))
    obs, vis = np.asarray(obs), np.asarray(vis)
    n_cams, n_pts = len(poses), len(pts0)
    n_grid = 2 * N_CTRL_B**2
    base = 6 * n_cams + 3 * n_pts + 1

    # Ideal (pre-distortion) normalized coords from the reference geometry —
    # the domain for the reference-model residual and the map comparison.
    ideal = np.empty_like(obs)
    for ci, (rv, tv) in enumerate(poses):
        m = vis[:, 0] == ci
        q = Rotation.from_rotvec(rv).apply(pts0[vis[m, 1]]) + tv
        ideal[m] = q[:, :2] / q[:, 2:3]
    ref_res = np.linalg.norm(ref_map(ideal) - obs, axis=1)
    obs_r = np.maximum(np.abs(obs[:, 0]) / (W / 2), np.abs(obs[:, 1]) / (H / 2))
    print(
        f"subset: {n_cams} images, {n_pts} points, {len(obs)} observations, "
        f"footprint to r={obs_r.max():.2f}; reference-model reproj rms "
        f"{np.sqrt((ref_res**2).mean()):.3f} px"
    )

    def unpack(x):
        p = x[: 6 * n_cams].reshape(n_cams, 6)
        pts = x[6 * n_cams : 6 * n_cams + 3 * n_pts].reshape(n_pts, 3)
        return p, pts, x[base - 1], x[base:]

    xy_lin = np.stack(
        np.meshgrid(np.linspace(-1, 1, 9), np.linspace(-1, 1, 9)), axis=-1
    ).reshape(-1, 2)
    des_lin = design_2d(xy_lin, N_CTRL_B, -1, 1)
    lin_basis = np.stack([np.ones(len(xy_lin)), xy_lin[:, 0], xy_lin[:, 1]], axis=1)
    proj_lin = np.linalg.lstsq(lin_basis, des_lin, rcond=None)[0]

    def residuals(x):
        p, pts, f, grid = unpack(x)
        res = np.empty((len(obs), 2))
        for ci in range(n_cams):
            m = vis[:, 0] == ci
            if m.any():
                res[m] = project(p[ci, :3], p[ci, 3:], pts[vis[m, 1]], f, grid) - obs[m]
        lin = np.concatenate(
            [proj_lin @ grid[: N_CTRL_B**2], proj_lin @ grid[N_CTRL_B**2 :]]
        )
        return np.concatenate([res.ravel(), 50.0 * lin])

    x0 = np.concatenate(
        [
            np.concatenate([np.concatenate(p) for p in poses]),
            pts0.ravel(),
            [f_ref * 1.02],
            np.zeros(n_grid),
        ]
    )
    n_par = len(x0)
    rows, cols = [], []
    for oi, (ci, pi) in enumerate(vis):
        cs = (
            list(range(6 * ci, 6 * ci + 6))
            + list(range(6 * n_cams + 3 * pi, 6 * n_cams + 3 * pi + 3))
            + list(range(base - 1, n_par))
        )
        for r in (2 * oi, 2 * oi + 1):
            rows += [r] * len(cs)
            cols += cs
    n_res = 2 * len(obs) + 6
    for r in range(2 * len(obs), n_res):
        rows += [r] * n_grid
        cols += list(range(base, n_par))
    spar = scipy.sparse.coo_matrix(
        (np.ones(len(rows)), (rows, cols)), shape=(n_res, n_par)
    )
    sol = scipy.optimize.least_squares(
        residuals, x0, jac_sparsity=spar, tr_solver="lsmr", xtol=1e-10, max_nfev=60
    )
    p, pts, f, grid = unpack(sol.x)

    # Structure drift vs the reference reconstruction (similarity-aligned).
    mu_a, mu_b = pts.mean(0), pts0.mean(0)
    a, b = pts - mu_a, pts0 - mu_b
    s = np.sqrt((b**2).sum() / (a**2).sum())
    u_, _, vt = np.linalg.svd(a.T @ b)
    rot = (u_ @ vt).T
    struct = np.sqrt(((s * a @ rot.T - b) ** 2).sum(axis=1).mean())

    sub = ideal[:: max(1, len(ideal) // 4000)]
    gt_map = ref_map(sub)
    u_lin = sub * f
    got = u_lin + grid_displace(np.clip(u_lin / np.array([W / 2, H / 2]), -1, 1), grid)
    merr, scale = map_error(got, gt_map)
    rms_reproj = np.sqrt((sol.fun[: 2 * len(obs)] ** 2).mean())
    print(
        f"grid solve: reproj {rms_reproj:.3f} px | f {f:.2f} (ref {f_ref:.2f}; "
        f"f carries the linear gauge and the scene-scale gauge, map scale "
        f"{scale:.4f}) | struct drift rmse {struct:.4f} | camera-map err vs "
        f"solved model rms {np.sqrt((merr**2).mean()):.3f} max {merr.max():.3f} px"
    )
    # Field scale for context: the reference model's displacement magnitude.
    ref_field = np.linalg.norm(gt_map - sub * f_ref, axis=1)
    print(
        f"reference distortion field over footprint: rms "
        f"{np.sqrt((ref_field**2).mean()):.3f} max {ref_field.max():.3f} px"
    )

    W, H = w_save, h_save


if __name__ == "__main__":
    experiment_a()
    experiment_b()
    experiment_c()
    experiment_d()
