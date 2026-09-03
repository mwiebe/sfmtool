# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""The depth reading: a position is stated only where one was measured."""

import numpy as np
import pytest

from seed_relax import Options, depth, release, structure

WIDTH, HEIGHT = 800, 600
FOCAL = 500.0
NAMES = ["cam/000.jpg", "cam/001.jpg", "cam/002.jpg"]
#: The frames sit on a straight line, a unit apart, so every point is seen
#: over the same two-unit baseline and its parallax is its depth's alone.
CENTRES = np.stack([np.array([float(f), 0.0, 0.0]) for f in range(len(NAMES))])
#: The pixel error every observation carries, alternating in sign so it is a
#: residual and not a shift of the principal point.
JITTER_PX = 0.4


def _camera():
    from sfmtool._sfmtool.geometry import CameraIntrinsics

    return CameraIntrinsics.from_dict(
        {
            "model": "SIMPLE_PINHOLE",
            "width": WIDTH,
            "height": HEIGHT,
            "parameters": {
                "focal_length": FOCAL,
                "principal_point_x": WIDTH / 2.0,
                "principal_point_y": HEIGHT / 2.0,
            },
        }
    )


def _grid(n, half, depth_z):
    """``n`` points on a square lattice at one depth, in world coordinates."""
    side = int(np.ceil(np.sqrt(n)))
    k = np.arange(n)
    x = -half + 2.0 * half * (k % side) / max(1, side - 1)
    y = -half + 2.0 * half * (k // side) / max(1, side - 1)
    return np.stack([x, y, np.full(n, -float(depth_z))], axis=1)


#: A near cloud with real parallax (11 degrees over the baseline) and a far one
#: with almost none (0.3 degrees), both carrying the same pixel error.
NEAR = _grid(30, 2.0, 10.0)
FAR = _grid(9, 8.0, 400.0)
#: One far point whose observations are exact.  Its own median residual is
#: zero, so only the member's floor states an error for it at all.
LUCKY = np.array([[0.0, 20.0, -400.0]])
POINTS = np.concatenate([NEAR, FAR, LUCKY])
N_NEAR, N_FAR = len(NEAR), len(FAR)
I_LUCKY = len(POINTS) - 1


def _project(cam, points, centre):
    return np.asarray(
        cam.ray_to_pixel_batch(np.ascontiguousarray(points - centre)), float
    )


def _member(points=POINTS, lucky=I_LUCKY):
    """A member whose every observation is off by :data:`JITTER_PX`."""
    import seed_candidate_eval as EV

    cam = _camera()
    obs_c, obs_i, obs_f, uv = [], [], [], []
    for f in range(len(NAMES)):
        px = _project(cam, points, CENTRES[f])
        for c in range(len(points)):
            obs_c.append(c)
            obs_i.append(f)
            obs_f.append(f * 1000 + c)
            off = 0.0 if c == lucky else JITTER_PX * (1.0 if (c + f) % 2 else -1.0)
            uv.append(px[c] + np.array([off, 0.0]))
    return EV.Member(
        0,
        "rotation_only",
        list(NAMES),
        cam,
        FOCAL,
        np.zeros((len(NAMES), 3)),
        np.zeros((len(NAMES), 3)),
        np.ones(len(NAMES), bool),
        np.asarray(points, float),
        (
            np.array(obs_c, np.int64),
            np.array(obs_i, np.int64),
            np.array(uv, float),
            np.array(obs_f, np.int64),
        ),
    )


def _state(points=POINTS):
    return {
        "frames": np.arange(len(NAMES), dtype=np.int64),
        "clusters": np.arange(len(points), dtype=np.int64),
        "quats": np.tile(np.array([1.0, 0.0, 0.0, 0.0]), (len(NAMES), 1)),
        "trans": np.stack([-c for c in CENTRES]),
        "points": np.asarray(points, float).copy(),
        "at_inf": np.zeros(len(points), bool),
    }


@pytest.fixture(name="read")
def _read():
    m = _member()
    st = _state()
    out, census = depth.demote_uncertain(m, m.camera, st, FOCAL)
    return m, st, out, census


def _demoted(before, after):
    return np.asarray(after["at_inf"], bool) & ~np.asarray(before["at_inf"], bool)


# ---------------------------------------------------------------- the rule


def test_a_wide_parallax_point_keeps_its_depth(read):
    _m, st, out, _census = read
    dem = _demoted(st, out)
    assert not dem[:N_NEAR].any()
    assert np.array_equal(out["points"][:N_NEAR], st["points"][:N_NEAR])


def test_a_low_parallax_far_point_loses_it(read):
    _m, st, out, census = read
    dem = _demoted(st, out)
    assert dem[N_NEAR:].all()
    assert census["n_demoted"] == N_FAR + 1
    assert census["n_finite_before"] == len(POINTS)
    assert census["n_finite_after"] == N_NEAR
    assert census["n_bearings_after"] == N_FAR + 1


def test_the_error_is_read_at_the_points_own_depth(read):
    """What is compared is the implied depth error, not the relative one.

    The far cloud's relative uncertainty is under the bound stated at the
    confident population's depth, so a single scalar bound would keep it. It
    is the far cloud's OWN depth that carries the same relative error past the
    support radius, and that is what the rule reads."""
    _m, st, out, census = read
    assert census["u_p90"] < census["support_scalar"]
    assert _demoted(st, out)[N_NEAR:].all()


def test_the_residual_floor_is_the_members_own_median(read):
    m, st, out, census = read
    assert census["eps_member_px"] == pytest.approx(JITTER_PX)
    # The lucky point's own observations are exact, so without the floor it
    # would state no error at all; with it, it is demoted like its neighbours.
    rows, slot_i, slot_c = structure.state_rows(m, st)
    rot, cen = structure.centres_of(st)
    resid = structure.reprojection(
        m.camera,
        rot,
        cen,
        st["points"],
        st["at_inf"],
        m.obs_uv[rows],
        slot_i,
        slot_c,
    )
    assert resid[slot_c == I_LUCKY] == pytest.approx(0.0, abs=1e-9)
    assert _demoted(st, out)[I_LUCKY]


def test_the_bound_is_read_on_the_confident_half(read):
    m, st, _out, census = read
    rows, slot_i, slot_c = structure.state_rows(m, st)
    rot, cen = structure.centres_of(st)
    uv = m.obs_uv[rows]
    rays = structure.world_rays(m.camera, rot, uv, slot_i)
    order, counts, starts = groups = depth.groups(slot_c, len(POINTS))
    theta = depth.widest_ray_angles(rays, *groups)
    resid = structure.reprojection(
        m.camera, rot, cen, st["points"], st["at_inf"], uv, slot_i, slot_c
    )
    eps = np.fmax(
        depth.point_medians(resid, order, counts, starts), census["eps_member_px"]
    )
    u = depth.uncertainties(theta, eps, FOCAL)
    conf = u <= np.median(u)
    assert int(conf.sum()) == census["n_confident"]
    # Every confident point is a near one: the far cloud is what the reading
    # is uncertain about.
    assert conf[N_NEAR:].sum() == 0
    assert census["support_r12"] == pytest.approx(
        depth.support_radius(st["points"][conf])
    )
    assert census["n_confident"] < census["n_finite_read"]
    # The scalar the census also states is that same radius carried to the
    # depth the confident points sit at.
    # The depth is read along the ray, so a point off the axis sits a little
    # further than its own Z.
    z_conf = float(np.median(np.abs(st["points"][conf][:, 2])))
    assert census["depth_med_confident"] == pytest.approx(z_conf, rel=0.05)
    assert census["support_scalar"] == pytest.approx(
        census["support_r12"] / census["depth_med_confident"]
    )


def test_a_demoted_point_keeps_its_observations_and_becomes_a_bearing(read):
    m, st, out, _census = read
    dem = _demoted(st, out)
    # Nothing is deleted: the same points, the same rows, nothing pruned.
    assert len(out["points"]) == len(out["at_inf"]) == len(POINTS)
    assert not len(structure.pruned_rows(out))
    before, _bi, _bc = structure.state_rows(m, st)
    after, _ai, _ac = structure.state_rows(m, out)
    assert np.array_equal(before, after)
    # A demoted point is a unit direction, and it points where its rays do.
    assert np.allclose(np.linalg.norm(out["points"][dem], axis=1), 1.0)
    towards = st["points"][dem] - CENTRES.mean(axis=0)
    towards /= np.linalg.norm(towards, axis=1, keepdims=True)
    assert np.all(np.einsum("ij,ij->i", out["points"][dem], towards) > 0.999)
    # The state it was handed is not the state it hands back.
    assert not np.asarray(st["at_inf"], bool).any()


def test_the_census_states_the_demotion_by_angle_bin(read):
    _m, _st, _out, census = read
    bins = {r["bin"]: r for r in census["bins"]}
    assert sum(r["n"] for r in census["bins"]) == len(POINTS)
    assert sum(r["demoted"] for r in census["bins"]) == census["n_demoted"]
    # The far cloud sits under a degree and goes whole; the near one is over
    # eight degrees and stays.
    assert bins["under 1 deg"]["n"] == N_FAR + 1
    assert bins["under 1 deg"]["frac"] == 1.0
    assert bins["over 8 deg"]["n"] == N_NEAR
    assert bins["over 8 deg"]["demoted"] == 0


# ------------------------------------------------------------- the primitives


def test_the_widest_angle_is_the_widest_pair():
    rays = np.array(
        [
            [0.0, 0.0, -1.0],
            [np.sin(np.radians(10.0)), 0.0, -np.cos(np.radians(10.0))],
            [np.sin(np.radians(30.0)), 0.0, -np.cos(np.radians(30.0))],
            [0.0, 0.0, -1.0],
        ]
    )
    slot_c = np.array([0, 0, 0, 1], np.int64)
    theta = depth.widest_ray_angles(rays, *depth.groups(slot_c, 2))
    assert np.degrees(theta[0]) == pytest.approx(30.0)
    # A single ray subtends nothing, which is what the floor rule already read.
    assert theta[1] == 0.0


def test_a_points_reading_is_the_median_over_its_own_rows():
    slot_c = np.array([0, 0, 0, 1, 2], np.int64)
    v = np.array([1.0, 9.0, 5.0, 4.0, 2.0])
    got = depth.point_medians(v, *depth.groups(slot_c, 4))
    assert got[:3].tolist() == [5.0, 4.0, 2.0]
    assert np.isnan(got[3])


def test_the_uncertainty_is_the_pixel_error_over_the_parallax():
    theta = np.array([np.radians(30.0), 0.0])
    u = depth.uncertainties(theta, np.array([2.0, 2.0]), 100.0)
    assert u[0] == pytest.approx(2.0 / (100.0 * 0.5))
    assert np.isnan(u[1])


def test_the_support_radius_is_the_twelfth_neighbours_distance():
    line = np.stack([np.arange(40.0), np.zeros(40), np.zeros(40)], axis=1)
    # On a unit line the twelfth nearest neighbour of an interior point is six
    # steps away, and the median over the set is that.
    assert depth.support_radius(line) == pytest.approx(6.0)
    assert depth.support_radius(line, k=2) == pytest.approx(1.0)


def test_the_support_radius_is_unstated_below_its_own_neighbourhood():
    assert depth.support_radius(np.zeros((depth.K_SUPPORT, 3))) is None
    assert depth.support_radius(np.zeros((depth.K_SUPPORT + 1, 3))) == 0.0


def test_the_neighbourhood_size_carries_its_provenance():
    assert depth.K_SUPPORT == 12
    assert "k_neighbors=12" in depth.K_SUPPORT_PROVENANCE["source_arg"]
    assert depth.K_SUPPORT_PROVENANCE["rule"]


# ---------------------------------------------------------------- refusals


def test_a_state_with_no_finite_point_refuses():
    m = _member()
    st = _state()
    st["at_inf"] = np.ones(len(POINTS), bool)
    out, census = depth.demote_uncertain(m, m.camera, st, FOCAL)
    assert out is st
    assert census["refused"] == "the state states no finite point"


def test_a_member_without_an_equivalent_focal_refuses():
    m = _member()
    st = _state()
    out, census = depth.demote_uncertain(m, m.camera, st, float("nan"))
    assert out is st
    assert "equivalent focal" in census["refused"]


def test_a_confident_half_smaller_than_the_neighbourhood_refuses():
    pts = np.concatenate([NEAR[:6], FAR[:4]])
    m = _member(pts, lucky=-1)
    st = _state(pts)
    out, census = depth.demote_uncertain(m, m.camera, st, FOCAL)
    assert out is st
    assert census["refused"] == "the confident half holds fewer than 13 points"
    assert census["n_demoted"] == 0


# --------------------------------------------------------------- determinism


def test_the_stage_is_a_function_of_its_inputs():
    got = []
    for _ in range(2):
        m = _member()
        out, census = depth.demote_uncertain(m, m.camera, _state(), FOCAL)
        got.append((out, census))
    for key in ("points", "at_inf"):
        assert (
            np.asarray(got[0][0][key]).tobytes() == np.asarray(got[1][0][key]).tobytes()
        )
    assert got[0][1] == got[1][1]


def test_the_pair_batch_size_cannot_change_the_reading(monkeypatch):
    m = _member()
    want, census = depth.demote_uncertain(m, m.camera, _state(), FOCAL)
    for chunk in (1, 4, 16):
        monkeypatch.setattr(depth, "PAIR_CHUNK", chunk)
        got, got_census = depth.demote_uncertain(m, m.camera, _state(), FOCAL)
        assert got["points"].tobytes() == want["points"].tobytes()
        assert got["at_inf"].tobytes() == want["at_inf"].tobytes()
        assert got_census == census


# ------------------------------------------------------------ the chain seam

CHAIN_FRAMES = 6
#: The chain's own frames, a quarter unit apart on a line.
CHAIN_CENTRES = np.stack([np.array([0.25 * f, 0.0, 0.0]) for f in range(CHAIN_FRAMES)])
CHAIN_STRIDE = 1000
#: Bearings the rotation-only model explained, near points it refused, and a
#: middle band whose parallax clears the angular floor by a hair, which the
#: re-estimation calls finite and whose depth the reading then withdraws.
CHAIN_FAR = _grid(60, 3000.0, 5000.0)
CHAIN_NEAR = _grid(60, 3.0, 10.0)
CHAIN_MID = _grid(40, 40.0, 150.0)
CHAIN_POINTS = np.concatenate([CHAIN_FAR, CHAIN_NEAR, CHAIN_MID])


def _chain_member():
    """A rotation-only member whose every observation is off by a jitter."""
    import seed_candidate_eval as EV

    cam = _camera()
    obs_c, obs_i, obs_f, uv, keep = [], [], [], [], []
    for f in range(CHAIN_FRAMES):
        px = _project(cam, CHAIN_POINTS, CHAIN_CENTRES[f])
        ok = np.isfinite(px).all(axis=1)
        for c in np.nonzero(ok)[0]:
            c = int(c)
            obs_c.append(c)
            obs_i.append(f)
            obs_f.append(f * CHAIN_STRIDE + c)
            off = JITTER_PX * (1.0 if (c + f) % 2 else -1.0)
            uv.append(px[c] + np.array([off, 0.0]))
            # The model kept the bearings and refused everything nearer.
            keep.append(c < len(CHAIN_FAR))
    dirs = CHAIN_POINTS - CHAIN_CENTRES.mean(axis=0)
    dirs = dirs / np.linalg.norm(dirs, axis=1, keepdims=True)
    return EV.Member(
        0,
        "rotation_only",
        [f"cam/{f:03d}.jpg" for f in range(CHAIN_FRAMES)],
        cam,
        FOCAL,
        np.zeros((CHAIN_FRAMES, 3)),
        np.zeros((CHAIN_FRAMES, 3)),
        np.ones(CHAIN_FRAMES, bool),
        dirs,
        (
            np.array(obs_c, np.int64),
            np.array(obs_i, np.int64),
            np.array(uv, float),
            np.array(obs_f, np.int64),
        ),
        keep=np.array(keep, bool),
    )


def _chain_source():
    """The handle the member was drawn from: its own clusters and no others."""
    import types

    cam = _camera()
    starts, images, features, affines = [0], [], [], []
    for c in range(len(CHAIN_POINTS)):
        for f in range(CHAIN_FRAMES):
            px = _project(cam, CHAIN_POINTS[c : c + 1], CHAIN_CENTRES[f])[0]
            if not np.isfinite(px).all():
                continue
            images.append(f)
            features.append(f * CHAIN_STRIDE + c)
            affines.append(np.array([[2.0, 0.0, px[0]], [0.0, 2.0, px[1]]]))
        starts.append(len(images))
    return types.SimpleNamespace(
        image_names=[f"cam/{f:03d}.jpg" for f in range(CHAIN_FRAMES)],
        refine_radius=8.0,
        cluster_starts=np.array(starts, np.int64),
        member_images=np.array(images, np.int64),
        member_features=np.array(features, np.int64),
        member_affines=np.stack(affines),
    )


@pytest.fixture(name="chain", scope="module")
def _chain():
    from unittest.mock import patch

    from seed_relax import pipeline

    on = pipeline.run_member(_chain_member(), _chain_source(), Options())
    # The chain as the reading was handed it: the same run with the stage
    # stubbed out to hand its state straight back, so the comparison below is
    # against the state the reading acted on.
    with patch.object(
        depth,
        "depth_stage",
        lambda mx, cam, st, f_eq, trace=None: (st, {"refused": "stubbed"}),
    ):
        off = pipeline.run_member(_chain_member(), _chain_source(), Options())
    return on, off


def test_the_reading_closes_the_chain(chain):
    on, _off = chain
    assert on.ok
    d = on.census["depth"]
    assert d.get("refused") is None
    assert d["n_demoted"] > 0
    # The counts the record ships are the ones the reading left behind.
    assert on.census["n_finite_final"] == d["n_finite_after"]
    assert on.census["n_infinity_final"] == d["n_bearings_after"]
    assert on.census["n_points_final"] == len(on.state["at_inf"])


def test_the_middle_band_is_what_the_reading_withdraws(chain):
    on, off = chain
    clusters = np.asarray(off.state["clusters"], np.int64)
    finite_off = set(clusters[~np.asarray(off.state["at_inf"], bool)].tolist())
    clusters_on = np.asarray(on.state["clusters"], np.int64)
    finite_on = set(clusters_on[~np.asarray(on.state["at_inf"], bool)].tolist())
    near = set(range(len(CHAIN_FAR), len(CHAIN_FAR) + len(CHAIN_NEAR)))
    mid = set(range(len(CHAIN_FAR) + len(CHAIN_NEAR), len(CHAIN_POINTS)))
    # The re-estimation graduated both bands; the reading keeps the near one
    # and withdraws the middle one.
    assert len(mid & finite_off) > 0.5 * len(mid)
    assert len(mid & finite_on) == 0
    assert len(near & finite_on) > 0.9 * len(near)


def test_a_stage_that_states_no_reading_is_the_chain_from_before_it(chain):
    _on, off = chain
    assert off.census["depth"] == {"refused": "stubbed"}
    inf = np.asarray(off.state["at_inf"], bool)
    assert off.census["n_finite_final"] == int((~inf).sum())
    # A stage that stated no reading states nothing in the released file's own
    # metadata either.
    assert "depth" not in release.tool_options(off, 0)


def test_the_manifest_carries_the_readings_census(chain):
    on, off = chain
    block = release.relaxation_block(on)["depth"]
    assert block["n_demoted"] == on.census["depth"]["n_demoted"]
    assert block["k_support"] == depth.K_SUPPORT
    assert block["support_r12"] > 0.0
    assert len(block["bins"]) == len(depth.ANGLE_NAMES)
    assert release.relaxation_block(off)["depth"]["refused"] == "stubbed"
    # And the released file names what the stage did.
    assert release.tool_options(on, 0)["depth"].startswith(
        f"{on.census['depth']['n_demoted']} demoted"
    )
    assert release.tool_options(on, 0)["points_finite"] == str(
        on.census["depth"]["n_finite_after"]
    )


def test_a_withdrawn_point_still_reaches_the_writer_as_a_bearing(chain):
    on, _off = chain
    data = {
        "names": list(on.member.names),
        "dims": [(WIDTH, HEIGHT)] * len(on.member.names),
        "obs_c": on.member.obs_c,
        "obs_i": on.member.obs_i,
        "obs_f": on.member.obs_f,
        "obs_uv": on.member.obs_uv,
        "n_img": len(on.member.names),
        "n_cl": on.member.n_cl,
    }
    _data_x, pts, keep, _res, at_inf = release.relaxed_arrays(on, data)
    mid = np.arange(len(CHAIN_FAR) + len(CHAIN_NEAR), len(CHAIN_POINTS))
    written = np.unique(on.member.obs_c[keep])
    kept_mid = np.intersect1d(mid, written)
    assert len(kept_mid) > 0
    assert at_inf[kept_mid].all()
    assert np.allclose(np.linalg.norm(pts[kept_mid], axis=1), 1.0)
