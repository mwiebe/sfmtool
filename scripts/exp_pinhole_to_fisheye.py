# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment F: can a pinhole start evolve into a fisheye, center-out?

On a fisheye capture (Kerry Park left fisheye), admit cluster observations
in growing radial rings. Two modes per ring:

- `pinhole`: plain affine factorization, no correction — the control.
- `evolve`:  alternate factorization with an L0 radial-profile refit,
  warm-started from the previous ring (the center-out curriculum).

Per ring we report the factorization residual, the inlier fraction, and the
evolved profile's error against the rig-calibrated OPENCV_FISHEYE's
equidistant→perspective correction (mod linear). At the last ring we try
the model-class upgrade: fit `s·f·tan(r/f)` to the evolved `r + g(r)` and
compare the implied focal with the calibrated one.

The perspective chart caps the experiment at θ ≈ 85° (beyond 90° there is
no perspective radius at all) — that boundary is itself part of the answer.

Run: pixi run -e dev python scripts/exp_pinhole_to_fisheye.py e_kerry_ws
"""

import json
import sys
from pathlib import Path

import numpy as np
import scipy.interpolate as si

WS = Path(sys.argv[1] if len(sys.argv) > 1 else "e_kerry_ws")
MIN_SPAN = 5
RINGS = [0.25, 0.35, 0.45, 0.55, 0.65]
TAU_PX = 3.0  # inlier threshold for the "explained" fraction
RNG = np.random.default_rng(13)


def design_1d(x, n_ctrl, lo, hi, degree=3):
    t = np.concatenate(
        [[lo] * degree, np.linspace(lo, hi, n_ctrl - degree + 1), [hi] * degree]
    )
    x = np.clip(x, lo, hi - 1e-9 * (hi - lo))
    return si.BSpline.design_matrix(x, t, degree).toarray()


# ── Data loading (same shapes as exp_cluster_bootstrap) ─────────────────────


def load_observations():
    from sfmtool._sfmtool.io import read_matches, read_sift, read_sift_metadata

    patches = sorted((WS / "matches").glob("*-clusters-patches.matches"))
    data = read_matches(patches[0])
    names = list(data["image_names"])
    counts = data["feature_counts"]
    prefix = data["metadata"]["workspace"]["contents"]["feature_prefix_dir"]
    positions, dims = [], []
    for i, name in enumerate(names):
        rel = Path(name)
        sp = WS / rel.parent / prefix / f"{rel.name}.sift"
        meta = read_sift_metadata(sp)["metadata"]
        dims.append((meta["image_width"], meta["image_height"]))
        s = read_sift(sp)
        positions.append(np.ascontiguousarray(s["positions_xy"][: int(counts[i])]))

    starts = np.asarray(data["cluster_starts"])
    mi = np.asarray(data["member_images"])
    mf = np.asarray(data["member_features"])
    st = np.asarray(data["member_status"])
    half = np.array([[w / 2, h / 2] for (w, h) in dims])
    rmax = float(np.linalg.norm(half[0]))

    obs_c, obs_i, obs_u = [], [], []
    n_cl = 0
    for c in range(len(starts) - 1):
        lo, hi = starts[c], starts[c + 1]
        sel = np.nonzero((st[lo:hi] == 0) | (st[lo:hi] == 1))[0] + lo
        if len(sel) < MIN_SPAN:
            continue
        for k in sel:
            obs_c.append(n_cl)
            obs_i.append(mi[k])
            obs_u.append(positions[mi[k]][mf[k]] - half[mi[k]])
        n_cl += 1
    return (
        np.asarray(obs_c),
        np.asarray(obs_i),
        np.asarray(obs_u, dtype=np.float64),
        rmax,
        len(names),
        n_cl,
    )


# ── Factorization (as in exp_cluster_bootstrap, trimmed) ────────────────────


def als_factorize(obs_c, obs_i, u_corr, n_img, n_cl, rounds=25, trim=0.05):
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
            sol = np.linalg.lstsq(xh, u_corr[s], rcond=None)[0]
            m_cam[i] = sol[:3].T
            t_cam[i] = sol[3]
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
    return res, keep


# ── Radial profile on a per-ring domain ──────────────────────────────────────


class Profile:
    """g(r)·rhat correction with knots on [0, domain]; constant beyond."""

    def __init__(self, knots=None, domain=1.0):
        self.knots = knots
        self.domain = domain

    def g(self, r):
        if self.knots is None:
            return np.zeros_like(r)
        rc = np.minimum(r, self.domain)
        des = design_1d(rc / self.domain, len(self.knots), 0.0, 1.0)
        return des @ self.knots

    def correct(self, u):
        r = np.linalg.norm(u, axis=1)
        rhat = np.divide(u, np.maximum(r, 1e-9)[:, None])
        return u + self.g(r)[:, None] * rhat


def fit_profile(u, target_delta, domain, n_knots):
    r = np.linalg.norm(u, axis=1)
    rhat = np.divide(u, np.maximum(r, 1e-9)[:, None])
    g_target = np.einsum("ij,ij->i", target_delta, rhat)
    des = design_1d(r / domain, n_knots, 0.0, 1.0)
    knots = np.linalg.lstsq(des, g_target, rcond=None)[0]
    gvals = des @ knots
    alpha = (gvals @ r) / (r @ r)  # remove the r-linear (gauge) part
    knots -= np.linalg.lstsq(des, alpha * r, rcond=None)[0]
    return Profile(knots, domain)


# ── Reference (rig calibration) ──────────────────────────────────────────────


def reference_correction(r_px):
    """Calibrated OPENCV_FISHEYE's equidistant→perspective correction (px)."""
    ref = json.loads((WS / "ref_camera.json").read_text())["parameters"]
    f = ref["fx"]
    ks = np.array([ref["k1"], ref["k2"], ref["k3"], ref["k4"]])
    thetad = r_px / f
    theta = thetad.copy()
    for _ in range(30):
        t2 = theta * theta
        poly = 1 + t2 * (ks[0] + t2 * (ks[1] + t2 * (ks[2] + t2 * ks[3])))
        dpoly = (
            3 * ks[0] * t2 + 5 * ks[1] * t2**2 + 7 * ks[2] * t2**3 + 9 * ks[3] * t2**4
        )
        theta -= (theta * poly - thetad) / (poly + dpoly)
    return f * np.tan(np.minimum(theta, 1.48)) - r_px, f


def mod_linear(g, r):
    return g - (g @ r) / (r @ r) * r


# ── Main ─────────────────────────────────────────────────────────────────────


def main():
    obs_c, obs_i, obs_u, rmax, n_img, n_cl = load_observations()
    obs_r = np.linalg.norm(obs_u, axis=1)
    print(
        f"{WS}: {n_img} images, {n_cl} clusters (span >= {MIN_SPAN}), "
        f"{len(obs_c)} observations, rmax {rmax:.0f} px"
    )
    print(
        f"{'mode':8s} {'ring':>5s} {'obs':>6s} {'resid':>7s} {'inlier<3px':>10s} "
        f"{'profile err (scale)':>22s}"
    )

    for mode in ["pinhole", "evolve"]:
        profile = Profile()
        for rho in RINGS:
            sel = obs_r < rho * rmax
            u = obs_u[sel]
            c, i = obs_c[sel], obs_i[sel]
            # Reindex clusters compactly for the factorization.
            uniq, c2 = np.unique(c, return_inverse=True)
            u_corr = profile.correct(u) if mode == "evolve" else u
            res, keep = als_factorize(c2, i, u_corr, n_img, len(uniq))
            if mode == "evolve":
                for _ in range(2):  # inner alternation at this ring
                    target = profile.correct(u)[keep] + res[keep] - u[keep]
                    profile = fit_profile(
                        u[keep], target, rho * rmax, 4 + int(round(8 * rho))
                    )
                    u_corr = profile.correct(u)
                    res, keep = als_factorize(c2, i, u_corr, n_img, len(uniq))
            rn = np.linalg.norm(res, axis=1)
            inlier = float((rn < TAU_PX).mean())

            r_eval = np.linspace(0.1, 0.98, 50) * rho * rmax
            g_ref, f_ref = reference_correction(r_eval)
            g_ref = mod_linear(g_ref, r_eval)
            if mode == "evolve":
                g_est = mod_linear(profile.g(r_eval), r_eval)
            else:
                g_est = np.zeros_like(r_eval)
            perr = np.sqrt(((g_est - g_ref) ** 2).mean())
            pscale = np.sqrt((g_ref**2).mean())
            print(
                f"{mode:8s} {rho:5.2f} {sel.sum():6d} "
                f"{np.sqrt((rn[keep] ** 2).mean()):7.2f} {inlier:10.3f} "
                f"{perr:9.2f} ({pscale:6.2f}) px"
            )

        if mode == "evolve":
            # Model-class upgrade: fit s·f·tan(r/f) to r + g(r) over the last
            # ring's domain; report the implied fisheye focal.
            r_fit = np.linspace(0.1, 0.95, 80) * RINGS[-1] * rmax
            mapped = r_fit + profile.g(r_fit)
            best = (np.inf, None, None)
            for f in np.linspace(60, 400, 341):
                basis = f * np.tan(np.clip(r_fit / f, 0, 1.45))
                s = (basis @ mapped) / (basis @ basis)
                err = np.sqrt(((s * basis - mapped) ** 2).mean())
                if err < best[0]:
                    best = (err, f, s)
            err, f_est, s = best
            # Pinhole is the f→inf limit; compare against a plain linear fit.
            s_lin = (r_fit @ mapped) / (r_fit @ r_fit)
            err_lin = np.sqrt(((s_lin * r_fit - mapped) ** 2).mean())
            print(
                f"model upgrade: tan-fit f = {f_est:.1f} px (calibrated fx "
                f"{f_ref:.1f}), fit err {err:.2f} px vs pinhole(linear) err "
                f"{err_lin:.2f} px"
            )


if __name__ == "__main__":
    main()
