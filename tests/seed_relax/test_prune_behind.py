# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""The cheirality refusal read per observation, and what the writer then sees.

One wrong match in a revisit frame used to demote a whole track to a bearing.
The estimate now drops that observation instead, and the state records it, so
neither an adjustment nor the release reads it again.
"""

import types

import numpy as np
import pytest

from seed_relax import release, structure

F = 500.0
WIDTH, HEIGHT = 640, 480
#: Four frames on a short arc that all see the points, and a fifth REVISIT
#: frame sitting beyond them along its own viewing axis.
CENTRES = np.array(
    [
        [0.0, 0.0, 0.0],
        [0.25, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.75, 0.0, 0.0],
        [0.0, 0.0, -20.0],
    ]
)
REVISIT = 4
#: Cluster 0 is the specimen: the revisit frame's observation of it is the
#: optical axis, which the point sits behind.  The other two are seen by the
#: agreeing frames alone.
POINTS = np.array([[0.0, 0.0, -10.0], [1.0, -0.5, -9.0], [-1.2, 0.6, -11.0]])
SPECIMEN = 0
FLOOR = np.radians(0.05)


@pytest.fixture(name="cam")
def _cam():
    from sfmtool._sfmtool.geometry import CameraIntrinsics

    return CameraIntrinsics.from_dict(
        {
            "model": "SIMPLE_PINHOLE",
            "width": WIDTH,
            "height": HEIGHT,
            "parameters": {
                "focal_length": F,
                "principal_point_x": WIDTH / 2.0,
                "principal_point_y": HEIGHT / 2.0,
            },
        }
    )


def _rows(cam):
    """``(obs_c, obs_i, obs_f, obs_uv)`` of the synthetic admission."""
    obs_c, obs_i, obs_f, uv = [], [], [], []
    for f in range(len(CENTRES) - 1):
        px = np.asarray(
            cam.ray_to_pixel_batch(np.ascontiguousarray(POINTS - CENTRES[f])), float
        )
        for c in range(len(POINTS)):
            obs_c.append(c)
            obs_i.append(f)
            obs_f.append(f * 100 + c)
            uv.append(px[c])
    # The wrong match: the revisit frame states the specimen at its principal
    # point, a ray the point sits behind.
    obs_c.append(SPECIMEN)
    obs_i.append(REVISIT)
    obs_f.append(REVISIT * 100 + SPECIMEN)
    uv.append([WIDTH / 2.0, HEIGHT / 2.0])
    return (
        np.array(obs_c, np.int64),
        np.array(obs_i, np.int64),
        np.array(obs_f, np.int64),
        np.array(uv, float),
    )


@pytest.fixture(name="member")
def _member(cam):
    import seed_candidate_eval as EV

    obs_c, obs_i, obs_f, uv = _rows(cam)
    return EV.Member(
        0,
        "rotation_only",
        [f"cam/{f:03d}.jpg" for f in range(len(CENTRES))],
        cam,
        F,
        np.zeros((len(CENTRES), 3)),
        np.zeros((len(CENTRES), 3)),
        np.ones(len(CENTRES), bool),
        POINTS / np.linalg.norm(POINTS, axis=1, keepdims=True),
        (obs_c, obs_i, uv, obs_f),
        keep=np.ones(len(obs_c), bool),
    )


def _state():
    """The placed state: identity rotations at the arc's own centres."""
    return {
        "frames": np.arange(len(CENTRES), dtype=np.int64),
        "clusters": np.arange(len(POINTS), dtype=np.int64),
        "quats": np.tile(np.array([1.0, 0.0, 0.0, 0.0]), (len(CENTRES), 1)),
        "trans": -CENTRES.copy(),
        "points": POINTS.copy(),
        "at_inf": np.zeros(len(POINTS), bool),
    }


# ── The recorded rows ─────────────────────────────────────────────────────


def test_a_state_that_pruned_nothing_reads_as_pruning_nothing():
    st = _state()
    assert structure.pruned_rows(st).tolist() == []
    rows = np.arange(9)
    assert structure.drop_pruned(st, rows) is rows


def test_pruned_rows_accumulate_and_stay_unique():
    st = structure.with_pruned(_state(), np.array([7, 3]))
    assert structure.pruned_rows(st).tolist() == [3, 7]
    st = structure.with_pruned(st, np.array([3, 11]))
    assert structure.pruned_rows(st).tolist() == [3, 7, 11]
    # The state it came from is untouched.
    assert structure.drop_pruned(_state(), np.arange(12)).tolist() == list(range(12))


def test_the_admission_rows_leave_out_what_was_pruned(cam, member):
    st = _state()
    rows, slot_i, slot_c = structure.state_rows(member, st)
    gone = rows[slot_i == REVISIT]
    assert len(gone) == 1
    st = structure.with_pruned(st, gone)
    kept, kept_i, kept_c = structure.state_rows(member, st)
    assert len(kept) == len(rows) - 1
    assert int(gone[0]) not in kept.tolist()
    assert (kept_i != REVISIT).all()
    assert len(kept_i) == len(kept) and len(kept_c) == len(kept)


# ── The estimate ──────────────────────────────────────────────────────────


def test_the_wrong_observation_is_dropped_and_the_point_stays_finite(cam, member):
    st = _state()
    rows, slot_i, slot_c = structure.state_rows(member, st)
    uv = member.obs_uv[rows]
    args = (cam, st["quats"], st["trans"], uv, slot_i, slot_c, len(POINTS), FLOOR)

    # As it stands, the one wrong observation demotes the whole track.
    _pts, at_inf, census, pruned = structure.estimate_points(*args)
    assert at_inf[SPECIMEN]
    assert census["n_behind"] == 1
    assert not pruned.any()
    assert census["n_pruned_obs"] == 0

    pts, at_inf, census, pruned = structure.estimate_points(*args, prune_behind=True)
    assert not at_inf.any()
    assert np.allclose(pts, POINTS, atol=1e-9)
    assert census["n_behind"] == 0
    assert census["n_finite"] == len(POINTS)
    assert census["n_finite_pruned"] == 1
    assert census["n_pruned_obs"] == 1
    # It is the revisit frame's row that goes.
    assert pruned.sum() == 1
    assert slot_i[pruned][0] == REVISIT
    assert slot_c[pruned][0] == SPECIMEN


# ── What the writer sees ──────────────────────────────────────────────────


def _data(member):
    return {
        "names": list(member.names),
        "dims": [(WIDTH, HEIGHT)] * len(member.names),
        "obs_c": member.obs_c,
        "obs_i": member.obs_i,
        "obs_f": member.obs_f,
        "obs_uv": member.obs_uv,
        "adm_rank": np.arange(member.n_cl),
        "cl_quality": np.zeros(member.n_cl),
        "n_img": len(member.names),
        "n_cl": member.n_cl,
    }


def _released(cam, member, state):
    result = types.SimpleNamespace(member=member, state=state, camera=cam)
    return release.relaxed_arrays(result, _data(member))


def _wrong_row(member, state):
    """The member row carrying the revisit frame's wrong match."""
    rows, slot_i, _slot_c = structure.state_rows(member, state)
    return int(rows[slot_i == REVISIT][0])


def _demoted(state):
    """The state the whole-track refusal leaves: the specimen a bearing."""
    st = dict(state)
    pts = np.asarray(st["points"], float).copy()
    at_inf = np.asarray(st["at_inf"], bool).copy()
    d = pts[SPECIMEN] / np.linalg.norm(pts[SPECIMEN])
    pts[SPECIMEN] = d
    at_inf[SPECIMEN] = True
    st["points"], st["at_inf"] = pts, at_inf
    return st


def test_the_bearing_the_whole_track_rule_leaves_carries_the_wrong_match(cam, member):
    # A bearing projects as a direction, so the observation the estimate could
    # not explain still reprojects and the writer keeps it.
    st = _demoted(_state())
    wrong = _wrong_row(member, st)
    _out, _pts, keep, res, at_inf = _released(cam, member, st)
    assert at_inf[SPECIMEN]
    assert keep[wrong]
    assert np.isfinite(res[wrong])


def test_the_release_does_not_carry_the_observation_the_estimate_refused(cam, member):
    st = _state()
    wrong = _wrong_row(member, st)
    _out, _pts, keep_before, _res, _inf = _released(cam, member, st)

    st = structure.with_pruned(st, [wrong])
    out, pts, keep, res, at_inf = _released(cam, member, st)
    assert not keep[wrong]
    assert not np.isfinite(res[wrong])
    # Every other row is untouched, and the point it observed is written as the
    # finite point its remaining observations state.
    others = np.ones(len(keep), bool)
    others[wrong] = False
    np.testing.assert_array_equal(keep[others], keep_before[others])
    assert not at_inf[SPECIMEN]
    assert np.allclose(pts[SPECIMEN], POINTS[SPECIMEN], atol=1e-9)
    assert len(out["obs_c"]) == len(member.obs_c)


def test_the_written_point_keeps_the_observations_it_was_solved_on(cam, member):
    st = _state()
    st = structure.with_pruned(st, [_wrong_row(member, st)])
    _out, pts, keep, _res, _inf = _released(cam, member, st)
    assert SPECIMEN in release.alive_clusters(pts, keep, member.obs_c).tolist()
    assert int(np.asarray(keep)[member.obs_c == SPECIMEN].sum()) == 4
