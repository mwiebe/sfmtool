# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""The reconciliation: several tracked points resting on one measurement."""

import types

import numpy as np
import pytest

from seed_relax import Options, pipeline, reconcile, structure

F = 500.0
WIDTH, HEIGHT = 640, 480
REFINE_RADIUS = 8.0
#: The member's own consensus tolerance: three pixels over its focal.
TOL_RAD = 3.0 / F
NAMES = ["cam/000.jpg", "cam/001.jpg", "cam/002.jpg"]
CENTRES = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
#: The world point every clean track in this module measures.
P = np.array([0.0, 0.0, -20.0])
#: Far enough that a one-unit baseline subtends less than the tolerance.
P_FAR = np.array([0.0, 0.0, -20000.0])


def _cam(focal=F):
    from sfmtool._sfmtool.geometry import CameraIntrinsics

    return CameraIntrinsics.from_dict(
        {
            "model": "SIMPLE_PINHOLE",
            "width": WIDTH,
            "height": HEIGHT,
            "parameters": {
                "focal_length": float(focal),
                "principal_point_x": WIDTH / 2.0,
                "principal_point_y": HEIGHT / 2.0,
            },
        }
    )


def _project(cam, point, image):
    xc = np.asarray(point, float)[None, :] - CENTRES[image][None, :]
    return np.asarray(cam.ray_to_pixel_batch(np.ascontiguousarray(xc)), float)[0]


def _build(rows, n_cl, points=None, at_inf=None):
    """``(member, state)`` over ``rows`` of ``(image, cluster, uv, scale)``."""
    import seed_candidate_eval as EV

    cam = _cam()
    obs_i = np.array([r[0] for r in rows], np.int64)
    obs_c = np.array([r[1] for r in rows], np.int64)
    uv = np.array([r[2] for r in rows], float)
    obs_f = np.arange(len(rows), dtype=np.int64)
    shapes = np.stack([np.array([[r[3], 0.0], [0.0, r[3]]], float) for r in rows])
    pts = np.tile(P, (n_cl, 1)) if points is None else np.asarray(points, float)
    m = EV.Member(
        0,
        "rotation_only",
        NAMES,
        cam,
        F,
        np.zeros((len(NAMES), 3)),
        np.zeros((len(NAMES), 3)),
        np.ones(len(NAMES), bool),
        pts,
        (obs_c, obs_i, uv, obs_f),
        shapes=shapes,
    )
    state = {
        "frames": np.arange(len(NAMES), dtype=np.int64),
        "clusters": np.arange(n_cl, dtype=np.int64),
        "quats": np.tile(np.array([1.0, 0.0, 0.0, 0.0]), (len(NAMES), 1)),
        "trans": np.ascontiguousarray(-CENTRES),
        "points": np.asarray(pts, float).copy(),
        "at_inf": np.zeros(n_cl, bool) if at_inf is None else np.asarray(at_inf, bool),
    }
    return m, state


def _run(m, state):
    return reconcile.reconcile_points(m, m.camera, state, REFINE_RADIUS, TOL_RAD, F)


# ------------------------------------------------------------ the relation


def test_bit_identical_positions_are_one_detection():
    slot_i = np.array([0, 0, 1], np.int64)
    uv = np.array([[10.0, 20.0], [10.0, 20.0], [10.0, 20.0]])
    gid = reconcile.detection_groups(slot_i, uv)
    # One detection in image 0 carrying two rows, another in image 1.
    assert gid[0] == gid[1]
    assert gid[2] != gid[0]


def test_exact_duplicates_join_without_a_tolerance():
    slot_i = np.zeros(2, np.int64)
    slot_c = np.array([0, 1], np.int64)
    uv = np.array([[10.0, 20.0], [10.0, 20.0]])
    gid = reconcile.detection_groups(slot_i, uv)
    # The near relation returns nothing (the pair is not two distinct places),
    # and the tangle exists anyway.
    a, b = reconcile.near_rows(slot_i, slot_c, uv, np.array([8.0, 80.0]), REFINE_RADIUS)
    assert len(a) == 0
    tangle, n_tangles, shared = reconcile.tangle_components(slot_c, gid, a, b, 2)
    assert n_tangles == 1
    assert tangle.tolist() == [0, 0]
    assert shared.all()


def test_the_tolerance_is_a_fraction_of_the_smaller_unit_scale():
    slot_i = np.zeros(2, np.int64)
    slot_c = np.array([0, 1], np.int64)
    uv = np.array([[0.0, 0.0], [1.0, 0.0]])
    # Unit scale 2 px: the bar is 0.6 * 2 = 1.2 px, and 1 px is inside it.
    r = np.array([16.0, 16.0])
    assert len(reconcile.near_rows(slot_i, slot_c, uv, r, REFINE_RADIUS)[0]) == 1
    # The same separation between features an octave finer: 0.6 * 1 = 0.6 px,
    # and the same 1 px is now two different places.
    r = np.array([8.0, 8.0])
    assert len(reconcile.near_rows(slot_i, slot_c, uv, r, REFINE_RADIUS)[0]) == 0


def test_the_smaller_of_the_two_sets_the_bar():
    slot_i = np.zeros(2, np.int64)
    slot_c = np.array([0, 1], np.int64)
    uv = np.array([[0.0, 0.0], [1.0, 0.0]])
    # A coarse row's own bar would admit this, but the radius test refuses the
    # pair before the position test is reached.
    r = np.array([80.0, 8.0])
    assert len(reconcile.near_rows(slot_i, slot_c, uv, r, REFINE_RADIUS)[0]) == 0
    # With the radius condition satisfied it is the SMALLER unit scale that
    # decides: two 8 px rows are 0.6 px apart at most, two 16 px rows 1.2.
    r = np.array([8.0, 8.4])
    assert len(reconcile.near_rows(slot_i, slot_c, uv, r, REFINE_RADIUS)[0]) == 0
    r = np.array([16.0, 16.8])
    assert len(reconcile.near_rows(slot_i, slot_c, uv, r, REFINE_RADIUS)[0]) == 1


def test_the_radius_agreement_is_a_tenth_of_the_larger():
    slot_i = np.zeros(2, np.int64)
    slot_c = np.array([0, 1], np.int64)
    uv = np.array([[0.0, 0.0], [0.2, 0.0]])
    assert (
        len(
            reconcile.near_rows(
                slot_i, slot_c, uv, np.array([16.0, 17.5]), REFINE_RADIUS
            )[0]
        )
        == 1
    )
    assert (
        len(
            reconcile.near_rows(
                slot_i, slot_c, uv, np.array([16.0, 18.0]), REFINE_RADIUS
            )[0]
        )
        == 0
    )


def test_a_row_never_pairs_with_its_own_cluster():
    slot_i = np.zeros(2, np.int64)
    slot_c = np.zeros(2, np.int64)
    uv = np.array([[0.0, 0.0], [0.2, 0.0]])
    a, _b = reconcile.near_rows(
        slot_i, slot_c, uv, np.array([16.0, 16.0]), REFINE_RADIUS
    )
    assert len(a) == 0


def test_a_tangle_is_the_connected_component_over_the_relation():
    # Three points in one chain: A shares with B, B shares with C, and A and C
    # share nothing.
    slot_i = np.array([0, 0, 1, 1], np.int64)
    slot_c = np.array([0, 1, 1, 2], np.int64)
    uv = np.array([[5.0, 5.0], [5.0, 5.0], [9.0, 9.0], [9.0, 9.0]])
    gid = reconcile.detection_groups(slot_i, uv)
    tangle, n_tangles, _shared = reconcile.tangle_components(
        slot_c, gid, np.zeros(0, np.int64), np.zeros(0, np.int64), 4
    )
    assert n_tangles == 1
    assert tangle.tolist() == [0, 0, 0, -1]


def test_the_tangle_numbering_does_not_depend_on_the_pair_order():
    slot_i = np.array([0, 0, 1, 1], np.int64)
    slot_c = np.array([2, 3, 0, 1], np.int64)
    uv = np.array([[5.0, 5.0], [5.0, 5.0], [9.0, 9.0], [9.0, 9.0]])
    gid = reconcile.detection_groups(slot_i, uv)
    z = np.zeros(0, np.int64)
    tangle, n, _s = reconcile.tangle_components(slot_c, gid, z, z, 4)
    assert n == 2
    # Numbered in point order: the tangle of points 0 and 1 comes first.
    assert tangle.tolist() == [0, 0, 1, 1]


# ------------------------------------------------------- the classifications


def _shared_pair_rows(bad_image=None, bad_offset=0.0):
    """Two clusters sharing image 0's detection, each with its own second view.

    ``bad_offset`` moves cluster 1's own row off the truth, which is what makes
    the union contradict."""
    cam = _cam()
    p0 = _project(cam, P, 0)
    rows = [
        (0, 0, tuple(p0), 2.0),
        (1, 0, tuple(_project(cam, P, 1)), 2.0),
        (0, 1, tuple(p0), 2.0),
        (2, 1, tuple(_project(cam, P, 2) + np.array([bad_offset, 0.0])), 2.0),
    ]
    if bad_image is not None:
        rows.append(
            (
                bad_image,
                1,
                tuple(_project(cam, P, bad_image) + np.array([200.0, 0.0])),
                2.0,
            )
        )
    return rows


def test_a_clean_tangle_merges_into_one_point():
    m, st = _build(_shared_pair_rows(), 2)
    m2, st2, census = _run(m, st)
    assert census["n_tangles"] == 1
    assert (census["merged"], census["culled"], census["refused_tangles"]) == (1, 0, 0)
    # One point, three rays, and the duplicate row on the shared detection is
    # gone from the admission.
    assert np.asarray(st2["clusters"]).tolist() == [0]
    assert census["rows_deduped"] == 1
    rows, _si, slot_c = structure.state_rows(m2, st2)
    assert len(rows) == 3
    assert set(slot_c.tolist()) == {0}
    assert np.allclose(np.asarray(st2["points"])[0], P, atol=1e-6)


def test_a_thin_union_merges_as_one_bearing():
    cam = _cam()
    p0 = _project(cam, P_FAR, 0)
    rows = [
        (0, 0, tuple(p0), 2.0),
        (1, 0, tuple(_project(cam, P_FAR, 1)), 2.0),
        (0, 1, tuple(p0), 2.0),
        (2, 1, tuple(_project(cam, P_FAR, 2)), 2.0),
    ]
    m, st = _build(rows, 2, points=np.tile(P_FAR, (2, 1)))
    _m2, st2, census = _run(m, st)
    assert census["merged_bearing"] == 1
    assert census["merged"] == 0
    assert np.asarray(st2["at_inf"]).tolist() == [True]


def test_a_single_removal_that_reconciles_culls_that_member():
    m, st = _build(_shared_pair_rows(bad_image=1), 2)
    m2, st2, census = _run(m, st)
    assert (census["merged"], census["culled"], census["refused_tangles"]) == (0, 1, 0)
    # Cluster 1 is what the arithmetic voted against, and every row of it has
    # left the admission.
    assert np.asarray(st2["clusters"]).tolist() == [0]
    rows, _si, _sc = structure.state_rows(m2, st2)
    assert len(rows) == 2
    assert np.allclose(np.asarray(st2["points"])[0], P, atol=1e-6)


def test_a_tangle_no_removal_reconciles_is_left_alone():
    cam = _cam()
    p0 = _project(cam, P, 0)
    rows = [(0, c, tuple(p0), 2.0) for c in range(3)]
    for c, off in enumerate((150.0, -150.0, 90.0)):
        rows.append((1, c, tuple(_project(cam, P, 1) + np.array([off, 0.0])), 2.0))
        rows.append((2, c, tuple(_project(cam, P, 2) + np.array([off, 0.0])), 2.0))
    m, st = _build(rows, 3)
    m2, st2, census = _run(m, st)
    assert census["n_tangles"] == 1
    assert (census["merged"], census["merged_bearing"], census["culled"]) == (0, 0, 0)
    assert census["refused_tangles"] == 1
    # Nothing moved: the same points, the same rows, the same member.
    assert m2 is m
    assert st2 is st
    assert np.asarray(st2["clusters"]).tolist() == [0, 1, 2]


def test_the_dedup_keeps_the_row_of_the_largest_cluster():
    cam = _cam()
    p0 = _project(cam, P, 0)
    # Cluster 0 carries three rows, cluster 1 two: the shared detection in
    # image 0 is stated by cluster 0's row.
    rows = [
        (0, 0, tuple(p0), 2.0),
        (1, 0, tuple(_project(cam, P, 1)), 2.0),
        (2, 0, tuple(_project(cam, P, 2)), 2.0),
        (0, 1, tuple(p0), 2.0),
        (1, 1, tuple(_project(cam, P, 1)), 2.0),
    ]
    m, st = _build(rows, 2)
    m2, st2, census = _run(m, st)
    assert census["merged"] == 1
    gone = structure.pruned_rows(st2)
    # Two rows are the same measurement twice; the ones retired are cluster 1's.
    assert census["rows_deduped"] == 2
    assert np.asarray(m.obs_c)[gone].tolist() == [1, 1]
    rows_left, _si, _sc = structure.state_rows(m2, st2)
    assert sorted(rows_left.tolist()) == [0, 1, 2]


def test_a_merge_never_reads_one_measurement_twice():
    m, st = _build(_shared_pair_rows(), 2)
    m2, st2, _census = _run(m, st)
    rows, slot_i, _sc = structure.state_rows(m2, st2)
    uv = np.asarray(m2.obs_uv, float)[rows]
    keys = {(int(i), float(p[0]), float(p[1])) for i, p in zip(slot_i, uv)}
    assert len(keys) == len(rows)


# ------------------------------------------------------------- the refusals


@pytest.mark.parametrize(
    "kw,reason",
    [
        ({"shapes": None}, "the member states no affine shapes"),
        ({"refine_radius": 0.0}, "the source file states no refine radius"),
    ],
)
def test_the_stage_refuses_where_the_unit_scale_is_unstated(kw, reason):
    m, st = _build(_shared_pair_rows(), 2)
    if "shapes" in kw:
        m.obs_shape = None
        rr = REFINE_RADIUS
    else:
        rr = kw["refine_radius"]
    m2, st2, census = reconcile.reconcile_points(m, m.camera, st, rr, TOL_RAD, F)
    assert census["refused"] == reason
    assert m2 is m and st2 is st


def test_a_state_with_no_tangle_is_handed_straight_back():
    cam = _cam()
    rows = [
        (i, c, tuple(_project(cam, P + np.array([3.0 * c, 0.0, 0.0]), i)), 2.0)
        for c in range(2)
        for i in range(3)
    ]
    m, st = _build(rows, 2, points=np.stack([P, P + np.array([3.0, 0.0, 0.0])]))
    m2, st2, census = _run(m, st)
    assert census["n_tangles"] == 0
    assert census["merged"] == 0
    assert m2 is m and st2 is st


# ------------------------------------------------------------- determinism


def test_the_stage_is_a_function_of_its_inputs():
    got = []
    for _ in range(2):
        m, st = _build(_shared_pair_rows(bad_image=1), 2)
        m2, st2, census = _run(m, st)
        got.append((m2, st2, census))
    a, b = got[0], got[1]
    assert a[2] == b[2]
    assert np.asarray(a[0].obs_c).tobytes() == np.asarray(b[0].obs_c).tobytes()
    for key in ("clusters", "points", "at_inf", "pruned_rows"):
        assert np.asarray(a[1][key]).tobytes() == np.asarray(b[1][key]).tobytes()


# ------------------------------------------------------------- the kill switch


def test_the_kill_switch_reads_the_environment():
    assert reconcile.reconcile_on({}) is True
    assert reconcile.reconcile_on({"SFMTOOL_RELAX_RECONCILE": "1"}) is True
    assert reconcile.reconcile_on({"SFMTOOL_RELAX_RECONCILE": "0"}) is False
    assert reconcile.reconcile_on({"SFMTOOL_RELAX_RECONCILE": " 0 "}) is False


# ---------------------------------------------------------------- end to end

N_FRAMES = 6
E2E_CENTRES = np.stack([np.array([0.25 * f, 0.0, 0.0]) for f in range(N_FRAMES)])
FEATURE_STRIDE = 1000


def _spread(n, depth, seed):
    k = np.arange(n)
    x = (-0.6 + 1.2 * ((k * 7 + seed) % n) / max(1, n - 1)) * depth
    y = (-0.4 + 0.8 * ((k * 3 + seed) % n) / max(1, n - 1)) * depth
    return np.stack([x, y, np.full(n, -float(depth))], axis=1)


FAR = _spread(60, 5000.0, 1)
NEAR = _spread(60, 10.0, 2)
#: The clusters the admission never held, and, at the end of the block, the
#: ORIENTATION DUPLICATES: five of them re-state a cluster the member already
#: holds, at its own pixel and its own scale, which is what the stage exists
#: for.
EXTRA = _spread(30, 14.0, 3)
N_DUP = 5


def _e2e_project(cam, points, centre):
    return np.asarray(
        cam.ray_to_pixel_batch(np.ascontiguousarray(points - centre)), float
    )


def _e2e_member():
    import seed_candidate_eval as EV

    cam = _cam()
    pts = np.concatenate([FAR, NEAR])
    obs_c, obs_i, obs_f, uv, keep, shapes = [], [], [], [], [], []
    for f in range(N_FRAMES):
        px = _e2e_project(cam, pts, E2E_CENTRES[f])
        for c in np.nonzero(np.isfinite(px).all(axis=1))[0]:
            obs_c.append(int(c))
            obs_i.append(f)
            obs_f.append(f * FEATURE_STRIDE + int(c))
            uv.append(px[c])
            keep.append(int(c) < len(FAR))
            # The affine the source file states for the same row, so the
            # member and the clusters the fill-in brings in read one radius.
            s = 2.0 if int(c) % 2 == 0 else 1.0
            shapes.append(np.array([[s, 0.0], [0.0, s]]))
    dirs = pts - E2E_CENTRES.mean(axis=0)
    dirs = dirs / np.linalg.norm(dirs, axis=1, keepdims=True)
    return EV.Member(
        0,
        "rotation_only",
        [f"cam/{f:03d}.jpg" for f in range(N_FRAMES)],
        cam,
        F,
        np.zeros((N_FRAMES, 3)),
        np.zeros((N_FRAMES, 3)),
        np.ones(N_FRAMES, bool),
        dirs,
        (
            np.array(obs_c, np.int64),
            np.array(obs_i, np.int64),
            np.array(uv, float),
            np.array(obs_f, np.int64),
        ),
        shapes=np.stack(shapes),
        keep=np.array(keep, bool),
    )


def _e2e_source():
    cam = _cam()
    starts, images, features, affines = [0], [], [], []
    blocks = [
        (np.concatenate([FAR, NEAR]), 0, 2.0, False),
        (EXTRA, len(FAR) + len(NEAR), 1.0, False),
        # The duplicates: the first five NEAR points again, under their own
        # feature indexes, which is how SIFT's orientation bins reach the file.
        (NEAR[:N_DUP], len(FAR) + len(NEAR) + len(EXTRA), 2.0, True),
    ]
    for points, base, scale, dup in blocks:
        for c in range(len(points)):
            for f in range(N_FRAMES):
                px = _e2e_project(cam, points[c : c + 1], E2E_CENTRES[f])[0]
                if not np.isfinite(px).all():
                    continue
                images.append(f)
                features.append(f * FEATURE_STRIDE + base + c)
                s = scale if dup or c % 2 == 0 else 0.5 * scale
                affines.append(np.array([[s, 0.0, px[0]], [0.0, s, px[1]]]))
            starts.append(len(images))
    return types.SimpleNamespace(
        image_names=[f"cam/{f:03d}.jpg" for f in range(N_FRAMES)],
        refine_radius=REFINE_RADIUS,
        cluster_starts=np.array(starts, np.int64),
        member_images=np.array(images, np.int64),
        member_features=np.array(features, np.int64),
        member_affines=np.stack(affines),
    )


@pytest.fixture(name="relaxed")
def _relaxed():
    return pipeline.run_member(_e2e_member(), _e2e_source(), Options())


def test_the_seam_reconciles_the_duplicates_the_fill_in_brought_in(relaxed):
    assert relaxed.ok
    rec = relaxed.census["reconcile"]
    assert rec.get("refused") is None
    assert rec["unit_frac"] == reconcile.NEAR_UNIT_FRAC
    # Every duplicate is one tangle of two points, and every one of them merges.
    assert rec["n_tangles"] == N_DUP
    assert rec["n_points_in_tangle"] == 2 * N_DUP
    assert rec["merged"] + rec["merged_bearing"] == N_DUP
    assert rec["culled"] == 0
    assert rec["refused_tangles"] == 0
    assert rec["n_points_after"] == rec["n_points"] - N_DUP


def test_the_manifest_block_states_what_the_stage_decided(relaxed):
    from seed_relax import release

    block = release.relaxation_block(relaxed)["reconcile"]
    for key in ("n_tangles", "merged", "culled", "refused_tangles", "rows_deduped"):
        assert block[key] == relaxed.census["reconcile"][key]
    opts = release.tool_options(relaxed, 0)
    assert "merged" in opts["reconcile"]


def test_the_kill_switch_leaves_the_chain_as_it_was(monkeypatch):
    off = pipeline.run_member(_e2e_member(), _e2e_source(), Options(reconcile=False))
    assert off.census["reconcile"] == {"held": "SFMTOOL_RELAX_RECONCILE"}
    from seed_relax import release

    assert release.relaxation_block(off)["reconcile"] == {
        "held": "SFMTOOL_RELAX_RECONCILE"
    }
    assert "reconcile" not in release.tool_options(off, 0)

    # The stage held is the stage absent: with the seam stubbed out to hand its
    # inputs straight back, the chain produces the same arrays bit for bit.
    monkeypatch.setattr(
        reconcile,
        "reconcile_stage",
        lambda mx, cam, st, rr, tol, f_eq, trace=None: (mx, st, {"held": "stub"}),
    )
    absent = pipeline.run_member(_e2e_member(), _e2e_source(), Options())
    for key in ("frames", "clusters", "quats", "trans", "points", "at_inf"):
        assert (
            np.asarray(off.state[key]).tobytes()
            == np.asarray(absent.state[key]).tobytes()
        )
    assert (
        np.asarray(off.member.obs_c).tobytes()
        == np.asarray(absent.member.obs_c).tobytes()
    )
    assert release.tool_options(off, 0) == release.tool_options(absent, 0)
