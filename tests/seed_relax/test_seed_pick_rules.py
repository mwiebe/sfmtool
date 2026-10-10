# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""The seed's qualification rule and the rotation core's resection retry floor.

Both live in `scripts/exp_fast_seed.py` (specs/core/geometry/seed-hypothesis-loop.md
§ "Rank" and specs/core/geometry/rotation-locked-resection.md § "Callers").
"""

import numpy as np
import pytest

import exp_fast_seed as FS
from sfmtool.geometry import CameraIntrinsics


def _res(**over):
    res = {"kept": FS.COMMIT_MIN_KEPT, "reach": 0.5, "spread": 0.2, "flags": []}
    res.update(over)
    return res


@pytest.mark.parametrize("reach", [0.04, 0.29, 0.59, 0.60, 0.61, 1.0])
def test_qualification_does_not_read_reach(reach):
    assert FS.qualifies(_res(reach=reach))


def test_qualification_reads_the_posed_count_and_the_scan_spread():
    assert not FS.qualifies(_res(kept=FS.COMMIT_MIN_KEPT - 1))
    assert not FS.qualifies(_res(spread=0.049))
    assert FS.qualifies(_res(spread=0.05))


@pytest.mark.parametrize(
    "flag", ["vote_divergence", "flat_scan", "edge_scan", "near_static_seed"]
)
def test_a_blocking_flag_disqualifies(flag):
    assert not FS.qualifies(_res(flags=[flag]))


def test_narrow_reach_and_low_consensus_do_not_block():
    assert FS.qualifies(_res(reach=0.05, flags=["narrow_reach", "low_consensus"]))


@pytest.mark.parametrize(
    "n, floor", [(12, 4), (40, 4), (41, 5), (48, 5), (62, 7), (71, 8), (100, 10)]
)
def test_retry_floor_is_a_tenth_of_the_offered_observations(n, floor):
    assert FS.resect_retry_floor(n) == floor


def test_retry_floor_never_falls_with_more_observations():
    floors = [FS.resect_retry_floor(n) for n in range(0, 300)]
    assert all(a <= b for a, b in zip(floors, floors[1:]))
    assert min(floors) == FS.RESECT_RETRY_MIN_COUNT


def _camera():
    return CameraIntrinsics(
        "PINHOLE",
        640,
        480,
        {
            "focal_length_x": 500.0,
            "focal_length_y": 500.0,
            "principal_point_x": 320.0,
            "principal_point_y": 240.0,
        },
    )


def _observations(cam, n_front, n_behind, t, seed=0):
    """``n_front`` exact observations in front of a camera at identity rotation
    and translation ``t``, plus ``n_behind`` points on the reflections of their
    rays through the camera centre: the sign-blind rows fit them exactly, and
    the in-front test removes them from every kept set."""
    rng = np.random.default_rng(seed)
    n = n_front + n_behind
    uv = np.column_stack([rng.uniform(40, 600, n), rng.uniform(40, 440, n)])
    rays = np.asarray(cam.pixel_to_ray_batch(np.ascontiguousarray(uv)))
    rays /= np.linalg.norm(rays, axis=1, keepdims=True)
    depth = rng.uniform(4.0, 9.0, n)
    sign = np.r_[np.ones(n_front), -np.ones(n_behind)]
    p_cam = rays * (depth * sign)[:, None]
    return p_cam - np.asarray(t)[None], uv


def test_retry_accepts_fewer_survivors_than_the_ordinary_floor_when_they_clear_its_share():
    cam = _camera()
    t = np.array([0.3, -0.2, 0.5])
    pts, uv = _observations(cam, 8, 22, t)
    out = FS.resect_locked(cam, [1.0, 0.0, 0.0, 0.0], pts, uv)
    assert out is not None
    assert int(np.asarray(out["inliers"]).sum()) == 8
    assert np.allclose(np.asarray(out["translation"]), t, atol=1e-6)


def test_retry_refuses_survivors_under_its_floor():
    cam = _camera()
    t = np.array([0.3, -0.2, 0.5])
    pts, uv = _observations(cam, 3, 27, t)
    assert FS.resect_locked(cam, [1.0, 0.0, 0.0, 0.0], pts, uv) is None
