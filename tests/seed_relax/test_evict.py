# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""The hand-over: a coarse observation a finer tracked feature covers."""

import types

import numpy as np
import pytest

from seed_relax import Options, evict, rings
from seed_relax.fleet_constants import RING_RATIO_P1

REFINE_RADIUS = 8.0
NAMES = ["cam/000.jpg", "cam/001.jpg", "cam/002.jpg"]

#: ``(cluster, uv, scale)`` per feature, repeated on every image.  An isotropic
#: affine of scale ``s`` reads as radius ``refine_radius * s`` and as keypoint
#: scale ``s``, so the drawn footprint is ``2.5 s`` and the refine footprint
#: ``8 s`` -- the two units the rule has to keep apart.
FEATURES = [
    # A coarse feature with a fine one well inside its DRAWN disk (2.5 * 4 =
    # 10 px): the hand-over the stage exists for.
    (0, (100.0, 100.0), 4.0),
    (1, (103.0, 100.0), 1.0),
    # A coarse feature with a fine one outside the drawn disk but inside the
    # refine disk (8 * 4 = 32 px): the refine grid would retire it, the drawn
    # footprint does not.
    (2, (300.0, 300.0), 4.0),
    (3, (315.0, 300.0), 1.0),
    # A same-scale neighbour, right on top: no band separates them, so nothing
    # is retired however close it sits.
    (4, (500.0, 500.0), 4.0),
    (5, (502.0, 500.0), 3.0),
]


def _rows(images=(0, 1, 2)):
    """``(uv, r_px, sigma, slot_i, slot_c)`` over ``images`` x FEATURES."""
    uv, sig, img, cl = [], [], [], []
    for i in images:
        for c, p, s in FEATURES:
            uv.append(p)
            sig.append(s)
            img.append(i)
            cl.append(c)
    sig = np.asarray(sig, float)
    return (
        np.asarray(uv, float),
        REFINE_RADIUS * sig,
        sig,
        np.asarray(img, np.int64),
        np.asarray(cl, np.int64),
    )


def _flag(**kw):
    uv, r, sig, i, c = _rows(**kw)
    flag, census = evict.covered_by_finer(uv, r, sig, i, c, 2.0)
    return flag, census, c


# --------------------------------------------------------------- containment


def test_the_footprint_the_rule_reads_is_the_drawn_one():
    flag, _census, cl = _flag()
    # Cluster 0's fine neighbour sits 3 px away, inside 2.5 * 4; cluster 2's
    # sits 15 px away, outside it but well inside 8 * 4.
    assert flag[cl == 0].all()
    assert not flag[cl == 2].any()


def test_the_fine_feature_is_never_the_one_retired():
    flag, _census, cl = _flag()
    for fine in (1, 3, 5):
        assert not flag[cl == fine].any()


def test_a_same_scale_neighbour_retires_nothing():
    flag, _census, cl = _flag()
    assert not flag[cl == 4].any()
    # And it is the SCALE test that spares it, not the distance: 2 px is well
    # inside the drawn disk.
    uv, r, sig, i, c = _rows()
    assert np.linalg.norm(uv[c == 4][0] - uv[c == 5][0]) < 2.5 * sig[c == 4][0]


def test_the_rule_never_lets_a_feature_cover_itself():
    # One image, one feature: the only row in the disk is the row itself.
    uv, r, sig, i, c = _rows(images=(0,))
    flag, census = evict.covered_by_finer(uv[:1], r[:1], sig[:1], i[:1], c[:1], 2.0)
    assert not flag.any()
    assert census["n_pairs_contained"] == 0


def test_a_band_finer_is_the_ratio_the_band_grid_is_cut_on():
    edges = rings.octave_edges(RING_RATIO_P1)
    assert evict.octave_ratio(edges) == 2.0
    # Every grid the constant can produce states the same ratio.
    assert evict.octave_ratio(rings.octave_edges(0.03)) == 2.0
    with pytest.raises(ValueError):
        evict.octave_ratio([float("inf")])


def test_the_ratio_bar_is_read_where_the_band_edge_is():
    uv = np.array([[0.0, 0.0], [1.0, 0.0]])
    sig = np.array([4.0, 2.0])
    i = np.zeros(2, np.int64)
    c = np.array([0, 1], np.int64)
    # Exactly one band apart: the radius ratio is 2, which the rule admits.
    flag, _ = evict.covered_by_finer(uv, REFINE_RADIUS * sig, sig, i, c, 2.0)
    assert flag.tolist() == [True, False]
    # A hair under one band, and nothing is retired.
    sig_near = np.array([4.0, 2.0000001])
    flag, _ = evict.covered_by_finer(uv, REFINE_RADIUS * sig_near, sig_near, i, c, 2.0)
    assert not flag.any()


def test_the_radius_is_the_fill_ins_own_reading():
    shapes = np.array([[[3.0, 0.0], [0.0, 5.0]], [[1.0, 0.0], [0.0, 1.0]]])
    got = evict.feature_radius(shapes, REFINE_RADIUS)
    assert got.tolist() == [0.5 * REFINE_RADIUS * 8.0, REFINE_RADIUS]


# ------------------------------------------------------------------ the sweep


def test_a_point_under_two_observations_goes_whole():
    slot_c = np.array([0, 0, 0, 1, 1, 2, 2], np.int64)
    keep = np.array([True, False, True, True, False, True, True], bool)
    ko, kp = evict.two_observation_sweep(slot_c, keep, 3)
    # Point 1 falls to one observation, so it and its survivor go.
    assert kp.tolist() == [True, False, True]
    assert ko.tolist() == [True, False, True, False, False, True, True]


def test_the_bands_state_what_went_in_and_what_came_out():
    ring = np.array([0, 0, 1, 1, 2], np.int64)
    slot_c = np.array([0, 0, 1, 1, 2], np.int64)
    keep_obs = np.array([False, False, True, True, True], bool)
    keep_pt = np.array([False, True, True], bool)
    got = evict.band_census(ring, keep_obs, slot_c, keep_pt)
    assert got == [
        {"ring": 0, "obs": 2, "obs_kept": 0, "points": 1, "points_kept": 0},
        {"ring": 1, "obs": 2, "obs_kept": 2, "points": 1, "points_kept": 1},
        {"ring": 2, "obs": 1, "obs_kept": 1, "points": 1, "points_kept": 1},
    ]


# ------------------------------------------------------------- the kill switch


@pytest.mark.parametrize(
    "env,want",
    [
        ({}, True),
        ({"SFMTOOL_RELAX_EVICT": "1"}, True),
        ({"SFMTOOL_RELAX_EVICT": ""}, True),
        ({"SFMTOOL_RELAX_EVICT": "0"}, False),
        ({"SFMTOOL_RELAX_EVICT": " 0 "}, False),
    ],
)
def test_the_kill_switch_reads_the_environment(env, want):
    assert evict.evict_on(env) is want


def test_the_switch_is_an_option_the_chain_carries():
    assert Options().evict is True
    assert Options(evict=False).evict is False


# ---------------------------------------------------------------- the refusals


def _source(prefix="features/sift-x"):
    return types.SimpleNamespace(
        metadata={"workspace": {"contents": {"feature_prefix_dir": prefix}}},
        refine_radius=REFINE_RADIUS,
    )


def _member_stub():
    return types.SimpleNamespace(
        names=list(NAMES),
        obs_i=np.zeros(1, np.int64),
        obs_f=np.zeros(1, np.int64),
        obs_c=np.zeros(1, np.int64),
        obs_uv=np.zeros((1, 2)),
        obs_shape=np.zeros((1, 2, 2)),
        rows_all=np.zeros(1, np.int64),
    )


def _state():
    return {
        "frames": np.zeros(1, np.int64),
        "clusters": np.zeros(1, np.int64),
        "quats": np.array([[1.0, 0.0, 0.0, 0.0]]),
        "trans": np.zeros((1, 3)),
        "points": np.zeros((1, 3)),
        "at_inf": np.ones(1, bool),
    }


def test_without_a_workspace_the_stage_refuses_and_changes_nothing():
    m, st = _member_stub(), _state()
    out, census = evict.evict_stage(m, None, st, _source(), None, REFINE_RADIUS, 1.0)
    assert out is st
    assert census["refused"] == "no workspace given"


def test_without_a_feature_directory_the_stage_refuses():
    m, st = _member_stub(), _state()
    out, census = evict.evict_stage(
        m, None, st, _source(prefix=""), "C:/nowhere", REFINE_RADIUS, 1.0
    )
    assert out is st
    assert "feature directory" in census["refused"]


def test_an_unreadable_sift_file_is_a_refusal_and_not_a_guess():
    m, st = _member_stub(), _state()
    out, census = evict.evict_stage(
        m, None, st, _source(), "C:/no/such/workspace", REFINE_RADIUS, 1.0
    )
    assert out is st
    assert census["refused"] and "no workspace given" not in census["refused"]


def test_a_scale_that_could_not_be_read_retires_nothing():
    m, st = _member_stub(), _state()
    out, census = evict.evict_covered(
        m, None, st, np.array([np.nan]), REFINE_RADIUS, 1.0
    )
    assert out is st
    assert census["refused"] == "a keypoint scale could not be read"


def test_a_member_without_a_refine_radius_refuses():
    m, st = _member_stub(), _state()
    out, census = evict.evict_covered(m, None, st, np.array([1.0]), None, 1.0)
    assert out is st
    assert "refine radius" in census["refused"]


# ---------------------------------------------------------------- the layout


def test_the_sift_file_is_the_one_beside_the_image():
    got = evict.sift_path("C:/ws", "features/sift-x", "images/left/frame_1.jpg")
    assert got.as_posix() == "C:/ws/images/left/features/sift-x/frame_1.jpg.sift"
    # A rig keeps one feature directory per camera, so two cameras of one
    # capture resolve to two directories.
    other = evict.sift_path("C:/ws", "features/sift-x", "images/right/frame_1.jpg")
    assert other.parent != got.parent


# --------------------------------------------------------------- determinism


def test_two_passes_produce_the_same_flag():
    a, _c, _cl = _flag()
    b, _c2, _cl2 = _flag()
    assert a.tobytes() == b.tobytes()


def test_the_pair_batch_size_cannot_change_the_flag(monkeypatch):
    want, census, _cl = _flag()
    for chunk in (1, 2, 3, 7):
        monkeypatch.setattr(evict, "_PAIR_CHUNK", chunk)
        got, got_census = evict.covered_by_finer(*_rows(), 2.0)
        assert got.tobytes() == want.tobytes()
        assert got_census == census


def test_the_order_the_rows_arrive_in_cannot_change_the_flag():
    uv, r, sig, i, c = _rows()
    want, _ = evict.covered_by_finer(uv, r, sig, i, c, 2.0)
    order = np.arange(len(uv))[::-1]
    got, _ = evict.covered_by_finer(
        uv[order], r[order], sig[order], i[order], c[order], 2.0
    )
    assert got[np.argsort(order)].tobytes() == want.tobytes()


# --------------------------------------------------------------- end to end

WIDTH, HEIGHT = 800, 600
#: ``(image, cluster, uv, scale)``.  Clusters 0-5 are FEATURES on every image;
#: 6 is a coarse cluster on two images, covered on one of them, so the sweep
#: takes it; 7 is the single-observation cover that takes it there and is then
#: swept itself.
E2E = (
    [(i, c, p, s) for i in (0, 1, 2) for c, p, s in FEATURES]
    + [(0, 6, (700.0, 100.0), 4.0), (1, 6, (700.0, 100.0), 4.0)]
    + [(0, 7, (702.0, 100.0), 1.0)]
)
N_CL_E2E = 8


def _e2e_member(rows=E2E, n_cl=N_CL_E2E):
    import seed_candidate_eval as EV
    from sfmtool._sfmtool.geometry import CameraIntrinsics

    cam = CameraIntrinsics.from_dict(
        {
            "model": "SIMPLE_PINHOLE",
            "width": WIDTH,
            "height": HEIGHT,
            "parameters": {
                "focal_length": 500.0,
                "principal_point_x": WIDTH / 2.0,
                "principal_point_y": HEIGHT / 2.0,
            },
        }
    )
    obs_i, obs_c, obs_f, uv, shapes = [], [], [], [], []
    for k, (i, c, p, s) in enumerate(rows):
        obs_i.append(i)
        obs_c.append(c)
        obs_f.append(k)
        uv.append(p)
        shapes.append(np.array([[s, 0.0], [0.0, s]]))
    pts = np.tile(np.array([0.0, 0.0, -20.0]), (n_cl, 1))
    return EV.Member(
        0,
        "rotation_only",
        NAMES,
        cam,
        500.0,
        np.zeros((len(NAMES), 3)),
        np.zeros((len(NAMES), 3)),
        np.ones(len(NAMES), bool),
        pts,
        (
            np.array(obs_c, np.int64),
            np.array(obs_i, np.int64),
            np.array(uv, float),
            np.array(obs_f, np.int64),
        ),
        shapes=np.stack(shapes),
    )


def _e2e_state(m, n_cl=N_CL_E2E):
    return {
        "frames": np.arange(len(NAMES), dtype=np.int64),
        "clusters": np.arange(n_cl, dtype=np.int64),
        "quats": np.tile(np.array([1.0, 0.0, 0.0, 0.0]), (len(NAMES), 1)),
        "trans": np.stack([np.array([-float(f), 0.0, 0.0]) for f in range(len(NAMES))]),
        "points": np.asarray(m.pts, float).copy(),
        "at_inf": np.zeros(n_cl, bool),
    }


@pytest.fixture(name="handed_over")
def _handed_over():
    from seed_relax import structure

    m = _e2e_member()
    st = _e2e_state(m)
    rows, _si, _sc = structure.state_rows(m, st)
    sigma = np.asarray(m.obs_shape, float)[rows][:, 0, 0]
    out, census = evict.evict_covered(m, m.camera, st, sigma, REFINE_RADIUS, None)
    return m, st, out, census


def test_the_stage_retires_the_covered_coarse_rows(handed_over):
    _m, _st, _out, census = handed_over
    # Cluster 0 on all three images, and cluster 6 on the one image its cover
    # is on.
    assert census["n_obs_covered"] == 4
    assert census["n_points_dropped_all_covered"] == 1
    assert census["n_points_dropped_by_two_obs"] == 2
    assert census["n_points_kept"] == 5
    assert census["n_obs_evicted"] == 6
    assert census["n_obs_kept"] == census["n_rows"] - 6


def test_the_state_comes_back_without_the_retired_evidence(handed_over):
    from seed_relax import structure

    m, st, out, _census = handed_over
    assert np.asarray(out["clusters"]).tolist() == [1, 2, 3, 4, 5]
    assert len(out["points"]) == len(out["at_inf"]) == 5
    gone = structure.pruned_rows(out)
    assert set(np.asarray(m.obs_c, np.int64)[gone].tolist()) == {0, 6, 7}
    # Nothing downstream reads them again.
    rows, _si, _sc = structure.state_rows(m, out)
    assert not set(rows.tolist()) & set(gone.tolist())
    # The original state is not mutated.
    assert len(st["clusters"]) == N_CL_E2E


def test_the_adjustment_runs_with_the_lens_held(handed_over):
    _m, _st, out, census = handed_over
    assert census["adjusted"] is True
    assert census["ba_n_obs"] == census["n_obs_kept"]
    assert len(out["quats"]) == len(NAMES)


def test_the_stage_is_a_function_of_its_inputs():
    from seed_relax import structure

    got = []
    for _ in range(2):
        m = _e2e_member()
        st = _e2e_state(m)
        rows, _si, _sc = structure.state_rows(m, st)
        sig = np.asarray(m.obs_shape, float)[rows][:, 0, 0]
        out, _c = evict.evict_covered(m, m.camera, st, sig, REFINE_RADIUS, None)
        got.append(out)
    for key in ("frames", "clusters", "quats", "trans", "points", "at_inf"):
        assert np.asarray(got[0][key]).tobytes() == np.asarray(got[1][key]).tobytes()


#: The six three-observation clusters alone: nothing the sweep can take, so a
#: rule that fires nowhere has to hand the state straight back.
E2E_FULL = [(i, c, p, s) for i in (0, 1, 2) for c, p, s in FEATURES]


def test_nothing_covered_leaves_the_state_untouched():
    from seed_relax import structure

    m = _e2e_member(E2E_FULL, len(FEATURES))
    st = _e2e_state(m, len(FEATURES))
    rows, _si, _sc = structure.state_rows(m, st)
    # Every drawn footprint below a pixel: no centre falls inside another, so
    # the rule fires nowhere.
    sigma = np.full(len(rows), 0.1)
    out, census = evict.evict_covered(m, m.camera, st, sigma, REFINE_RADIUS, None)
    assert out is st
    assert census["n_obs_covered"] == 0
    assert census["n_obs_evicted"] == 0
    assert census["adjusted"] is False
