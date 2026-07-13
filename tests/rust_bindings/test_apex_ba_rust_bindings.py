# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Tests for the apex-solver bundle-adjustment EVALUATION binding
(``sfmtool._sfmtool.geometry.bundle_adjust_apex``, branch ``apex-ba-eval``).

Evaluation scaffolding: these tests exist to validate the binding used by
``scripts/exp_pinhole_bootstrap.py`` while the apex-solver crate is under
evaluation; they leave the tree with the binding if the crate is rejected.

Covers: pose round-trips through the apex SE3 (the axis-angle exp-map
convention of scipy's ``Rotation.from_rotvec``, reimplemented here in numpy
because the test env has no scipy), recovery of a perturbed synthetic
pinhole scene with the focal fixed and free (both via the custom
single-focal factor — the built-in ProjectionFactor's zero-residual
cheirality convention collapses anchor-free BAs into the scene reflection;
see the binding's module docs), and repeat-run determinism (recorded as a
warning rather than a failure if the solver turns out nondeterministic).
"""

import warnings

import numpy as np
import numpy.testing as npt

from sfmtool._sfmtool.geometry import bundle_adjust_apex

N_CAMS = 20
N_PTS = 500
F_TRUE = 800.0


def _rotvec_to_matrix(rv):
    """Rodrigues exp map, (N, 3) -> (N, 3, 3); scipy from_rotvec convention."""
    rv = np.atleast_2d(rv)
    theta = np.linalg.norm(rv, axis=1)
    out = np.zeros((len(rv), 3, 3))
    for k, (v, th) in enumerate(zip(rv, theta)):
        kx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        if th < 1e-12:
            out[k] = np.eye(3) + kx
        else:
            out[k] = (
                np.eye(3)
                + np.sin(th) / th * kx
                + (1 - np.cos(th)) / (th * th) * (kx @ kx)
            )
    return out


def _rot_angle(ra, rb):
    """Rotation angle between matrix batches (N, 3, 3).

    Frobenius form ``||Ra - Rb||_F = 2 sqrt(2) sin(theta/2)``: linear near
    zero, unlike the trace/arccos form whose precision floor is sqrt(eps)."""
    d = np.linalg.norm((ra - rb).reshape(len(ra), -1), axis=1)
    return 2 * np.arcsin(np.clip(d / (2 * np.sqrt(2)), 0.0, 1.0))


def _look_at_cam_from_world(eye, target, up):
    """cam_from_world (R, t) with +z forward, x_cam = R X + t."""
    fwd = target - eye
    fwd = fwd / np.linalg.norm(fwd)
    right = np.cross(fwd, up)
    right = right / np.linalg.norm(right)
    down = np.cross(fwd, right)
    r = np.column_stack([right, down, fwd]).T
    return r, -r @ eye


def _matrix_to_rotvec(r):
    """Inverse of the exp map for a single matrix (angle in [0, pi])."""
    angle = np.arccos(np.clip((np.trace(r) - 1) / 2, -1.0, 1.0))
    if angle < 1e-12:
        return np.zeros(3)
    if angle > np.pi - 1e-3:
        # Near pi the antisymmetric part vanishes; use (R + I)/2 ~ a a^T.
        b = (r + np.eye(3)) / 2
        k = int(np.argmax(np.diag(b)))
        axis = b[:, k] / np.sqrt(b[k, k])
        return angle * axis / np.linalg.norm(axis)
    axis = np.array([r[2, 1] - r[1, 2], r[0, 2] - r[2, 0], r[1, 0] - r[0, 1]])
    return angle * axis / np.linalg.norm(axis)


def _synthetic_scene(seed=7):
    """Ring of cameras around a point cloud; every point seen everywhere."""
    rng = np.random.default_rng(seed)
    pts = rng.uniform(-2.0, 2.0, size=(N_PTS, 3))
    rvec = np.zeros((N_CAMS, 3))
    tvec = np.zeros((N_CAMS, 3))
    for i in range(N_CAMS):
        ang = 2 * np.pi * i / N_CAMS
        eye = np.array([10 * np.cos(ang), 2.0 * np.sin(3 * ang), 10 * np.sin(ang)])
        r, t = _look_at_cam_from_world(eye, np.zeros(3), np.array([0.0, 1.0, 0.0]))
        rvec[i] = _matrix_to_rotvec(r)
        tvec[i] = t
    obs_i = np.repeat(np.arange(N_CAMS, dtype=np.uint32), N_PTS)
    obs_p = np.tile(np.arange(N_PTS, dtype=np.uint32), N_CAMS)
    rots = _rotvec_to_matrix(rvec)
    xc = np.einsum("nij,nj->ni", rots[obs_i], pts[obs_p]) + tvec[obs_i]
    assert (xc[:, 2] > 1.0).all()
    obs_xy = F_TRUE * xc[:, :2] / xc[:, 2:3]
    return rvec, tvec, pts, obs_p, obs_i, obs_xy


def _perturbed(rng, rvec, tvec, pts):
    return (
        rvec + rng.normal(0, 0.01, rvec.shape),
        tvec + rng.normal(0, 0.05, tvec.shape),
        pts + rng.normal(0, 0.02, pts.shape),
    )


def _rms_reproj(res, obs_p, obs_i, obs_xy):
    rots = _rotvec_to_matrix(np.asarray(res["rvec"]))
    xc = (
        np.einsum("nij,nj->ni", rots[obs_i], np.asarray(res["points"])[obs_p])
        + np.asarray(res["tvec"])[obs_i]
    )
    proj = res["f"] * xc[:, :2] / xc[:, 2:3]
    return float(np.sqrt(((proj - obs_xy) ** 2).mean()))


def _solve(scene, init, opt_f, f0, max_iterations=60):
    rvec0, tvec0, pts0 = init
    _, _, _, obs_p, obs_i, obs_xy = scene
    return bundle_adjust_apex(
        obs_p,
        obs_i,
        obs_xy,
        rvec0,
        tvec0,
        pts0,
        f0,
        opt_f,
        1.0,
        max_iterations,
    )


def test_pose_round_trip_exact():
    """From the exact solution, poses must round-trip through apex SE3
    (conversion fidelity + rotation convention vs the rotvec exp map)."""
    scene = _synthetic_scene()
    rvec, tvec, pts = scene[:3]
    res = _solve(scene, (rvec, tvec, pts), opt_f=False, f0=F_TRUE, max_iterations=1)
    # Zero residual at the optimum: the accepted step is ~0, so any pose
    # movement is conversion error.  Rotations compared as rotations (the
    # returned rotvec is canonical, magnitude <= pi).
    rot_err = _rot_angle(
        _rotvec_to_matrix(np.asarray(res["rvec"])), _rotvec_to_matrix(rvec)
    )
    assert rot_err.max() < 1e-9
    npt.assert_allclose(np.asarray(res["tvec"]), tvec, atol=1e-9)
    npt.assert_allclose(np.asarray(res["points"]), pts, atol=1e-9)
    assert res["f"] == F_TRUE
    assert _rms_reproj(res, *scene[3:]) < 1e-9


def test_recovers_scene_fixed_f():
    """Built-in ProjectionFactor path (opt_f=False) recovers the scene."""
    scene = _synthetic_scene()
    rng = np.random.default_rng(11)
    res = _solve(scene, _perturbed(rng, *scene[:3]), opt_f=False, f0=F_TRUE)
    assert res["iterations"] <= 60
    assert _rms_reproj(res, *scene[3:]) < 1e-3


def test_recovers_scene_free_f():
    """Custom single-focal factor path (opt_f=True) recovers scene + focal."""
    scene = _synthetic_scene()
    rng = np.random.default_rng(13)
    res = _solve(scene, _perturbed(rng, *scene[:3]), opt_f=True, f0=F_TRUE * 1.02)
    assert _rms_reproj(res, *scene[3:]) < 1e-3
    assert abs(res["f"] - F_TRUE) / F_TRUE < 1e-3  # within 0.1%


def test_repeat_run_determinism():
    """Two identical runs should be bitwise identical; if not, record it
    prominently instead of failing (an evaluation finding, not a bug in
    our binding)."""
    scene = _synthetic_scene()
    rng_a = np.random.default_rng(17)
    init = _perturbed(rng_a, *scene[:3])
    res_a = _solve(scene, init, opt_f=True, f0=F_TRUE * 1.02)
    res_b = _solve(scene, init, opt_f=True, f0=F_TRUE * 1.02)
    bitwise = (
        res_a["f"] == res_b["f"]
        and res_a["final_cost"] == res_b["final_cost"]
        and res_a["iterations"] == res_b["iterations"]
        and np.array_equal(np.asarray(res_a["rvec"]), np.asarray(res_b["rvec"]))
        and np.array_equal(np.asarray(res_a["tvec"]), np.asarray(res_b["tvec"]))
        and np.array_equal(np.asarray(res_a["points"]), np.asarray(res_b["points"]))
    )
    if not bitwise:
        warnings.warn(
            "EVALUATION FINDING: apex-solver repeat runs are NOT bitwise "
            f"identical (final_cost {res_a['final_cost']!r} vs "
            f"{res_b['final_cost']!r}, f {res_a['f']!r} vs {res_b['f']!r})",
            stacklevel=1,
        )
        # Still expect agreement to solver-noise level.
        npt.assert_allclose(
            np.asarray(res_a["points"]), np.asarray(res_b["points"]), atol=1e-6
        )
