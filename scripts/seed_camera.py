# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Camera context, cluster loading and `.sfmr` writers of the seed pipeline.

Holds the per-run camera context -- base model, focal, radial-spline
coefficients -- that every camera the seed builds is made from, the
model-generic depth and field tests taken against it, the loader that turns a
`*-clusters-patches.matches` file into flat observation arrays, the
densification that triangulates every cluster at final poses, and the writers
that put the result on disk (the final save and the debug seed snapshots).

`scripts/exp_fast_seed.py` installs the context through `set_camera_context`
and drives these writers; `scripts/check_fisheye_seed_primitives.py` exercises
the same primitives against a synthetic capture.
"""

import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from sfmtool._sfmtool.geometry import (
    CameraIntrinsics,
    refine_absolute_pose as _refine_absolute_pose,
    reprojection_residuals as _reprojection_residuals,
)

WS = Path(sys.argv[1] if len(sys.argv) > 1 else "e_seoul_ws")
MIN_SPAN_BA = int(  # min distinct images for a cluster to become a point
    os.environ.get("SFMTOOL_MIN_SPAN_BA", "2")
)
F_GRID = [0.55, 0.7, 0.9, 1.2, 1.6]  # focal candidates, in units of max(w, h)
TRIM_PX = 4.0  # BA inter-round observation trim threshold


# Fronto-parallel prior damping the seed writer's own tilt solve: the
# out-of-plane slope b_z is the weakly observed, depth-like direction, so an
# undamped solve on a narrow-baseline window is noise.  Exposed for sweeping.
WRITER_FRONTO_LAM = float(os.environ.get("SFMTOOL_WRITER_FRONTO_LAM", "0.3"))


# The whole script works in the CANONICAL camera frame (-Z forward, +Y up):
# poses are canonical world->camera, 3D points are world points, observations
# are FULL (un-centered) pixel coordinates, and every projection goes through
# the native `CameraIntrinsics` batch functions.  The world frame is the
# COLMAP-world gauge inherited from the affine factorization (irrelevant to
# the reprojection residuals and absorbed by the eval's similarity alignment);
# only the writer rotates it by W to reach the .sfmr canonical world.
_CAM_WH = (
    None  # (w, h) of the shared pinhole; set in main() from the uniform image dims
)


# ── Fisheye seed camera context (scripts/notes-fisheye-seed.md, Phase 1) ─────
#
# The finalization's twin of ``exp_fast_seed``'s camera context: a per-run
# (model, focal) pair behind every camera this script builds.  Default
# SIMPLE_PINHOLE — the code path this script has always run, byte-identical.
#
# The context is INSTALLED BY THE CALLER, never inferred here: the confirmed
# equidistant verdict lives in stage 1, and ``exp_fast_seed`` hands it over
# through ``set_camera_context`` before ``finalize_seed_from_dict``.  A
# standalone bootstrap run therefore stays pinhole unless the caller says
# otherwise.  As in stage 1, the fisheye model is EQUIDISTANT_FISHEYE (the
# SIMPLE_PINHOLE analog for `theta = r/f`, closed form both ways with an
# analytic pixel Jacobian), and only the primitives that build their camera
# through ``make_cam`` become equidistant — the photometric embed's patch
# geometry is Phase 4.
_CAM_CONTEXT = {
    "model": "SIMPLE_PINHOLE",
    "focal": None,
    "bspline": None,
    "theta_max": None,
}


def set_camera_context(model, focal=None, bspline=None, theta_max=None):
    """Install the per-run camera context (see the block comment above).

    ``bspline`` / ``theta_max`` are the radial-spline coefficients and the END
    OF THEIR DOMAIN, ignored by every other model.  The domain end is measured
    in the base model's own radial coordinate — incidence angle under the
    equidistant base, normalized image-plane radius under the pinhole one (see
    ``_SPLINE_MODEL``) — and the keyword keeps its fisheye name over both.  The
    finalization's spline rung is the only caller that ever passes them, and it
    passes whichever coordinate its own base measures; stage 1 stays
    base-model."""
    _CAM_CONTEXT["model"] = model
    _CAM_CONTEXT["focal"] = None if focal is None else float(focal)
    _CAM_CONTEXT["bspline"] = (
        None if bspline is None else np.ascontiguousarray(bspline, dtype=np.float64)
    )
    _CAM_CONTEXT["theta_max"] = None if theta_max is None else float(theta_max)


def camera_context():
    """The active ``(model, focal)`` context as a plain dict (a copy)."""
    return dict(_CAM_CONTEXT)


# ── The radial-spline models, keyed off the base model ──────────────────────
#
# SFMTOOL_FISHEYE and SFMTOOL_PINHOLE are ONE model with two bases: a monotone
# cubic B-spline over the base's own radial coordinate ``d``, added to ``d``
# before the focal scales it to pixels
# (specs/formats/sfmtool-camera-models.md).  ``d`` is the incidence angle
# ``theta`` under the equidistant base and the normalized image-plane radius
# ``rho = tan(theta)`` under the pinhole one, and that is the ONLY difference:
# same parameter head, same coefficients, same gauge, same monotonicity
# invariant.  Everything this script does with a spline camera is therefore
# written in ``d``, and this table is the single point where the base model
# picks which coordinate that is.
_SPLINE_MODEL = {
    # base model -> (spline model, domain-end parameter, fisheye base?)
    "SIMPLE_PINHOLE": ("SFMTOOL_PINHOLE", "bspline_rho_max", False),
    "SFMTOOL_PINHOLE": ("SFMTOOL_PINHOLE", "bspline_rho_max", False),
    "EQUIDISTANT_FISHEYE": ("SFMTOOL_FISHEYE", "bspline_theta_max", True),
    "SIMPLE_RADIAL_FISHEYE": ("SFMTOOL_FISHEYE", "bspline_theta_max", True),
    "SFMTOOL_FISHEYE": ("SFMTOOL_FISHEYE", "bspline_theta_max", True),
}


def spline_model(model=None):
    """``(spline model, domain-end parameter, fisheye base?)`` of the
    radial-spline promotion of ``model`` — the installed context's model when
    omitted.  A model the table does not carry is read as a fisheye base, the
    one this script promoted before the pinhole spline existed."""
    m = _CAM_CONTEXT["model"] if model is None else model
    return _SPLINE_MODEL.get(m, _SPLINE_MODEL["EQUIDISTANT_FISHEYE"])


def _d_of_theta(theta, fisheye_base=None):
    """The base model's radial coordinate at incidence angle ``theta``: the
    angle itself under a fisheye base, ``tan(theta)`` under a pinhole one."""
    if fisheye_base is None:
        fisheye_base = spline_model()[2]
    return theta if fisheye_base else np.tan(theta)


def _theta_of_d(d, fisheye_base=None):
    """The incidence angle at radial coordinate ``d`` — inverse of
    :func:`_d_of_theta`."""
    if fisheye_base is None:
        fisheye_base = spline_model()[2]
    return d if fisheye_base else np.arctan(d)


def fisheye_stage1():
    """Whether a FISHEYE-BASE camera context is installed — the finalization's
    twin of ``exp_fast_seed.fisheye_stage1``, and the single test every Phase-4
    branch is gated on.  A fisheye base only ever arrives through
    ``set_camera_context``, which stage 1 calls on a CONFIRMED both-cells verdict
    (routing by default, unless ``SFMTOOL_FISHEYE_SEED=0`` refuses it), so no
    capture the arbitration did not confirm as fisheye can reach any fisheye
    branch.

    Read off ``_SPLINE_MODEL`` rather than as "not the pinhole default": the
    spline rung promotes a PINHOLE capture too, and SFMTOOL_PINHOLE is a
    non-default model whose base is still the pinhole one.  What the branches
    gated here ask — the ray-range depth measure, the imaged-cone field test,
    the writer's model-generic surfel arm — is a question about the base, not
    about whether a distortion rung ran.  Every model reachable before that
    promotion answers exactly as the old test did."""
    return spline_model()[2]


# Absolute plausibility floor on a released focal, as a multiple of
# max(w, h).  The pinhole value is the long-standing bound; the equidistant
# one is the low end of the focal-vote kernel's FOV-derived band
# (specs/core/geometry/focal-vote.md), because `theta = r/f` ties focal to field of
# view — a >180 deg capture's own focal sits BELOW the pinhole floor (kerry:
# f ~ 138 px against 0.3 x 480 = 144), which would reject every honest solve.
_FOCAL_FLOOR_MULT = {
    "SIMPLE_PINHOLE": 0.3,
    "EQUIDISTANT_FISHEYE": 0.075,
    # The spline rung's promotion of the same lens: same floor as its base.
    "SFMTOOL_FISHEYE": 0.075,
    "SFMTOOL_PINHOLE": 0.3,
}


def focal_floor():
    """The context's absolute focal plausibility floor, px."""
    return _FOCAL_FLOOR_MULT.get(_CAM_CONTEXT["model"], 0.075) * max(_CAM_WH)


def make_cam(f=None):
    """The context camera at focal ``f`` (the context focal when omitted).

    SIMPLE_PINHOLE by default (principal point at the image centre);
    EQUIDISTANT_FISHEYE under an installed fisheye context, or the matching
    radial-spline model once the finalization's spline rung has promoted it
    (same map plus the context's spline).  All of them share the same three
    base parameters — a promoted model adds its distortion — so
    this builds one dict.  The images share one size (see main()), so one
    camera serves every projection; ``ray_to_pixel_batch`` /
    ``pixel_to_ray_batch`` map canonical camera-space points <-> full
    pixels."""
    w, h = _CAM_WH
    if f is None:
        f = _CAM_CONTEXT["focal"]
    params = {
        "focal_length": float(f),
        "principal_point_x": w / 2.0,
        "principal_point_y": h / 2.0,
    }
    model = _CAM_CONTEXT["model"]
    if model in ("SFMTOOL_FISHEYE", "SFMTOOL_PINHOLE"):
        coeffs = np.asarray(_CAM_CONTEXT["bspline"], dtype=np.float64)
        params[spline_model(model)[1]] = float(_CAM_CONTEXT["theta_max"])
        params["bspline_coeff_count"] = float(len(coeffs))
        for i, c in enumerate(coeffs):
            params[f"bspline_c{i}"] = float(c)
    return CameraIntrinsics.from_dict(
        {
            "model": model,
            "width": int(w),
            "height": int(h),
            "parameters": params,
        }
    )


def make_cam_bspline(f, coeffs, d_max):
    """The context camera promoted to its radial-spline model at
    ``(f, coeffs)`` on the spline domain ``[0, d_max]`` — SFMTOOL_FISHEYE
    under an equidistant base, SFMTOOL_PINHOLE under a pinhole one, with
    ``d_max`` in whichever radial coordinate that base measures
    (see ``_SPLINE_MODEL``).

    An all-zero ``coeffs`` is the base model's own map — bit for bit,
    projection, inverse and pixel Jacobian alike (the model's zero-spline
    identity) — so the promotion itself moves nothing; the coefficients are
    the dimensionless radial correction the base map cannot
    express, with ``f`` staying the central scale under the model's
    center-anchored gauge."""
    w, h = _CAM_WH
    model, d_key, _ = spline_model()
    cc = np.asarray(coeffs, dtype=np.float64)
    params = {
        "focal_length": float(f),
        "principal_point_x": w / 2.0,
        "principal_point_y": h / 2.0,
        d_key: float(d_max),
        "bspline_coeff_count": float(len(cc)),
    }
    for i, c in enumerate(cc):
        params[f"bspline_c{i}"] = float(c)
    return CameraIntrinsics.from_dict(
        {
            "model": model,
            "width": int(w),
            "height": int(h),
            "parameters": params,
        }
    )


def _cam_depth(p_cam):
    """Model-aware "distance in front" of CANONICAL camera-frame points.

    ``-z`` for the perspective family (the canonical camera looks down -Z), the
    ray RANGE ``|p|`` under a fisheye context.  Past 90 degrees off axis ``-z``
    is negative, so every ordering, median and ratio taken over it inverts on
    exactly the periphery a >180 degree capture exists to image; the range is
    the distance those observations actually sit at, and the two agree on axis.
    Broadcasts over the last axis, so it takes an ``(n, 3)`` block or one row."""
    if fisheye_stage1():
        return np.linalg.norm(p_cam, axis=-1)
    return -np.asarray(p_cam)[..., 2]


def _in_field(cam, p_cam):
    """Whether CANONICAL camera-frame points lie inside the cone ``cam`` images.

    Perspective family: the half-space in front, ``-z > 0`` — equivalently
    ``theta < 90 deg``.  A fisheye's imaged cone is not that half-space, so the
    same statement is made model-generically: ``theta <= r_max / f``, with
    ``r_max`` the INSCRIBED image-circle radius (half the smaller image
    dimension), i.e. the largest off-axis angle the sensor carries.  For a
    pinhole the two coincide; for the equidistant map the pinhole form would
    discard the whole 90-degree-plus annulus, which is 18-37% of every fisheye
    entry's detected features."""
    p = np.asarray(p_cam, dtype=np.float64)
    if not fisheye_stage1():
        return -p[..., 2] > 0
    rng = np.linalg.norm(p, axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        theta = np.arccos(np.clip(-p[..., 2] / np.maximum(rng, 1e-300), -1.0, 1.0))
    return (rng > 0) & (theta <= _field_theta_max(cam))


def _field_d_max(cam):
    """The largest radial coordinate ``cam`` images: the base model's own
    coordinate at the inscribed image circle's rim — the incidence angle
    ``theta`` under a fisheye base, the normalized image-plane radius
    ``rho = tan(theta)`` under a pinhole one.

    Distortion-free the pixel radius is ``f * d`` either way, so the rim is
    ``r_max / f`` outright.  Once a distortion rung has promoted the camera
    (the spline rung's SFMTOOL_FISHEYE / SFMTOOL_PINHOLE, or a legacy
    SIMPLE_RADIAL_FISHEYE) the pixel radius is ``f * (d + delta(d))``, not
    ``f * d``, so the rim is read back through the model's own inverse rather
    than divided out."""
    r_max = 0.5 * min(cam.width, cam.height)
    if cam.model in ("SIMPLE_RADIAL_FISHEYE", "SFMTOOL_FISHEYE", "SFMTOOL_PINHOLE"):
        cx, cy = cam.principal_point
        ray = cam.pixel_to_ray(cx + r_max, cy)
        theta = float(np.arccos(np.clip(-ray[2], -1.0, 1.0)))
        return _d_of_theta(theta, spline_model(cam.model)[2])
    return r_max / float(cam.focal_lengths[0])


# The fisheye entry point: under a fisheye base the radial coordinate IS the
# incidence angle, which is the reading ``_in_field`` and its callers take.
_field_theta_max = _field_d_max


def _colmap_proj_jacobian(cam, xc_col, s_flip, rel_step=1e-5):
    """Per-row 2x3 pixel Jacobian ``d(u, v)/d x`` of ``cam`` at the COLMAP-frame
    (+Z forward) camera-space points ``xc_col``.

    The camera object projects CANONICAL (-Z forward) rays, and the two frames
    differ by the involution ``S = diag(1, -1, -1)``, so ``x_can = S x_col`` and
    the chain rule is a column-wise sign flip: ``J_col = J_can · S``.  ``J_can``
    is a central difference of ``ray_to_pixel_batch``, which is what makes this
    model-generic — the same measure ``WarpMap`` takes its warp Jacobian by,
    valid for the equidistant map at every theta including past 90 degrees,
    where no image-plane form exists.  The step is relative to each row's own
    range so it is scale-free (the projection is degree-0 homogeneous in the
    ray, so only the direction matters).

    Returns an ``(n, 2, 3)`` array."""
    xc_can = np.ascontiguousarray(xc_col * s_flip)
    h = rel_step * np.maximum(np.linalg.norm(xc_can, axis=1), 1e-12)
    out = np.zeros((len(xc_can), 2, 3))
    for j in range(3):
        step = np.zeros_like(xc_can)
        step[:, j] = h
        plus = np.asarray(
            cam.ray_to_pixel_batch(np.ascontiguousarray(xc_can + step)),
            dtype=np.float64,
        )
        minus = np.asarray(
            cam.ray_to_pixel_batch(np.ascontiguousarray(xc_can - step)),
            dtype=np.float64,
        )
        out[:, :, j] = (plus - minus) / (2.0 * h)[:, None] * s_flip[j]
    return out


def reproj_res_one(cam, rvec_i, tvec_i, x_pts, uv, invalid=1e6):
    """(proj − obs) pixel residuals of one image's gathered points under a
    single canonical world->camera pose (rvec/tvec), via the native
    ``reprojection_residuals`` kernel.  Behind-camera observations get
    ``invalid`` on their x component (never an inlier), mirroring the old
    ``max(z, 1e-6)`` clamp.  Returns an (N, 2) array."""
    q = Rotation.from_rotvec(rvec_i).as_quat()[[3, 0, 1, 2]][None, :]
    n = len(uv)
    return _reprojection_residuals(
        cam,
        q,
        np.ascontiguousarray(tvec_i, dtype=np.float64)[None, :],
        np.ascontiguousarray(x_pts, dtype=np.float64),
        np.ascontiguousarray(uv, dtype=np.float64),
        np.zeros(n, np.uint32),
        np.arange(n, dtype=np.uint32),
        invalid,
    )


# ── Data loading ─────────────────────────────────────────────────────────────


def relative_warps(shapes, obs_c, reference_members):
    """Reference->member warps ``W = S.S_ref^-1`` from stored absolute shapes.

    Since .matches format version 5 the affine's leading 2x2 is the member's
    ABSOLUTE affine shape ``S = W.S_ref`` -- the detector's canonical unit
    frame mapped onto that member's image pixels -- not the reference-relative
    warp.  Consumers that want image-space SIZE read ``S`` directly; consumers
    that want the warp between two views of one patch (the surfel writer's
    tilt solve) invert the cluster's own reference row ``S_ref`` to get it
    back.

    ``shapes`` is ``(K, 2, 2)``, ``obs_c`` the per-member cluster id, and
    ``reference_members`` the per-cluster global member index (``0xFFFFFFFF``
    where the cluster carries no reference).  Members of a reference-less
    cluster have no recoverable warp -- a derived file whose reference fell
    outside its restriction -- and come back as zeros, exactly like the
    all-zero rows of members that were never evaluated.
    """
    refs = np.asarray(reference_members)
    have = refs != np.iinfo(np.uint32).max
    out = np.zeros_like(shapes)
    rows = np.nonzero(have[obs_c])[0]
    if len(rows):
        s_ref = shapes[refs[obs_c[rows]].astype(np.int64)]
        out[rows] = shapes[rows] @ np.linalg.inv(s_ref)
    return out


def load_clusters(matches_data=None, preselected=False):
    """Patch clusters as flat observation arrays with refined positions.

    Everything geometric comes straight from the .matches file: image
    dimensions from the images section and member positions from the stored
    affines' last column (the absolute refined keypoint position).  The
    admission itself — drop unrefinable clusters, keep reference/kept
    members, restrict to selected images, span-filter — is the
    matches-format crate's ``select_clusters`` derivation; this function
    just reshapes the selected file into the flat observation arrays.

    ``matches_data`` (optional): an already-open ``MatchesFile`` handle to
    reuse.  Default None: read the workspace's clusters-patches file.

    ``preselected``: the handle IS the admission — a restriction stage already
    ran the selection and wrote it, so this load reshapes its arrays as they
    stand.  Re-selecting such a file would drop every cluster whose reference
    member fell outside the restriction (the derived file records the
    absent-reference sentinel, which the selection reads as unrefinable), so a
    preselected load never re-selects.  Its cluster ids are the file's own,
    which is what the stages upstream of it name.

    Restricting to a subset of images is NOT an option here: it is the
    restriction stage's job, and its artifact is what a ``preselected`` load
    reads.  Two independent restrictions of one source renumber independently,
    which is the coupling that stage exists to remove.
    """
    _t_load = time.perf_counter()
    if matches_data is not None:
        data = matches_data
    else:
        from sfmtool._sfmtool.io import MatchesFile

        override = os.environ.get("SFMTOOL_MATCHES")
        patches = (
            [Path(override)]
            if override
            else sorted(WS.glob("matches/*-clusters-patches.matches"))
        )
        print(f"matches file: {patches[0]}")
        data = MatchesFile(patches[0])
    names = list(data.image_names)
    dims = [(int(w), int(h)) for w, h in np.asarray(data.image_dims)]

    # File-level selection (native): reference/kept members, clusters spanning
    # >= MIN_SPAN_BA distinct images.  Cluster order (by source id) and member
    # order are preserved, so the observation stream below matches the
    # selection's CSR layout directly.
    sel = data if preselected else data.select_clusters(min_span=MIN_SPAN_BA)
    starts = np.asarray(sel.cluster_starts, dtype=np.int64)
    sizes = np.diff(starts)
    n_cl = len(sizes)
    aff = np.asarray(sel.member_affines)

    obs_i = np.asarray(sel.member_images, dtype=np.int64)
    if preselected:
        n2r = int((sizes == 2).sum())
        print(
            f"preselected admission over {len(names)} images: span-2 {n2r}, "
            f"span>=3 {n_cl - n2r} usable clusters"
        )

    # Admission order (best first) — used for both the cap and the tiers:
    # highest span first (ties broken by cluster id for determinism).  The
    # selected file keeps at most one reference/kept member per (cluster,
    # image), so each cluster's span IS its member count.
    #
    # No admission cap: growth and triangulation see every usable cluster
    # (a capped set can disconnect a chain-shaped capture — south-building
    # fragmented at 36/128).  The ordering instead selects which clusters'
    # observations enter the BAs (the top MAX_CLUSTERS by adm_rank).
    order = np.lexsort((np.arange(n_cl), -sizes))
    adm_rank = np.empty(n_cl, dtype=np.int64)
    adm_rank[order] = np.arange(n_cl)

    print(f"load_clusters: {time.perf_counter() - _t_load:.2f}s")
    obs_c = np.repeat(np.arange(n_cl, dtype=np.int64), sizes)
    # The affine's last column is the member's absolute refined keypoint
    # position; the 2x2 block is its ABSOLUTE affine shape (`S_ref | x_ref`
    # for the reference row).  The surfel writer wants the reference-relative
    # warp, so recover it through each cluster's reference row.
    out = {
        "names": names,
        "dims": dims,
        "obs_c": obs_c,
        "obs_i": obs_i,
        "obs_f": np.asarray(sel.member_features, dtype=np.int64),
        "obs_uv": np.ascontiguousarray(aff[:, :, 2], dtype=np.float64),
        "obs_warp": np.ascontiguousarray(
            relative_warps(aff[:, :, :2], obs_c, sel.reference_members),
            dtype=np.float64,
        ),
        "obs_ref": np.asarray(sel.member_status) == 0,
        "adm_rank": adm_rank,
        # Worst (max) finite warp-consistency residual over the selected
        # members — lower is better; clusters where no member entered the
        # consistency fit rank last (inf).
        "cl_quality": np.asarray(sel.cluster_worst_consistency(), dtype=np.float64),
        "refine_radius": data.refine_radius,
        "n_img": len(names),
        "n_cl": n_cl,
    }
    return out


def p3p_resect(uv, x_pts, f0, wh):
    """Minimal-sample absolute pose: RANSAC P3P over 2D-3D candidates.

    The trimmed-LS ``pose_refine`` needs a decent inlier fraction; a
    junk-match-dominated image (dino img 52: ~7-10% true 2D-3D pairs from a
    4x physical scale gap) defeats it, while minimal 3-point sampling finds
    the consensus routinely.  Uses the native Lambda Twist estimator
    (specs/core/geometry/absolute-pose.md); a tight 4 px threshold matches the
    bootstrap's TRIM_PX (a loose consensus is mostly junk on a
    wrong-match-heavy image and anchoring the verification BA on it drags
    the pose).  ``uv`` are full pixels.  Returns (rvec, tvec, inlier mask
    over the given obs) or None."""
    from sfmtool._sfmtool.geometry import estimate_absolute_pose

    ans = estimate_absolute_pose(
        np.ascontiguousarray(uv),
        np.ascontiguousarray(x_pts),
        camera=make_cam(f0),
        max_error_px=4.0,
        seed=0,
    )
    if ans is None:
        return None
    # The estimator already returns a canonical world-to-camera pose, which
    # is the frame the whole script works in — no flip.
    q = np.asarray(ans["quaternion_wxyz"])
    rv = Rotation.from_quat(q[[1, 2, 3, 0]]).as_rotvec()
    tv = np.asarray(ans["translation"], dtype=np.float64)
    return rv, tv, np.asarray(ans["inliers"], dtype=bool)


def pose_refine(uv, x_pts, rv0, tv0, f):
    """Pose-only resection of one image against known 3D points.

    Trimmed iterations (native ``refine_absolute_pose``): repeatedly refit L2
    on the best-fitting 60% of the observations, then a final refit on the
    < 3 px inliers.  A plain L2 warm-up is dragged by the junk observations'
    leverage, and a robust loss has near-zero gradient when every residual
    starts as a 100 px "outlier" — trimming from a decent init has neither
    problem.  Canonical world->camera pose in, canonical pose out."""
    q0 = Rotation.from_rotvec(rv0).as_quat()[[3, 0, 1, 2]]
    out = _refine_absolute_pose(
        make_cam(f),
        np.ascontiguousarray(uv, dtype=np.float64),
        np.ascontiguousarray(x_pts, dtype=np.float64),
        q0,
        np.ascontiguousarray(tv0, dtype=np.float64),
        5,  # trim rounds
        0.6,  # keep fraction
        3.0,  # final inlier px
    )
    q = np.asarray(out["quaternion_wxyz"])
    rv = Rotation.from_quat(q[[1, 2, 3, 0]]).as_rotvec()
    tv = np.asarray(out["translation"], dtype=np.float64)
    return rv, tv, float(out["inlier_fraction"])


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


# ── Perspective conversion + triangulation ───────────────────────────────────


def dense_structure(all_c, all_i, all_u, f, rvec, tvec, pts, posed, quiet=False):
    """The fine structure pass the coarse growth defers.

    Growth and the BAs run on the capped cluster subset (the scaling lever), so
    the structure they leave behind covers only that subset — a frame whose
    covisibility lives entirely in capped-out clusters ends up posed but with
    zero observations (a pose-only skeleton). The poses are final here, so
    triangulate EVERY cluster at those poses and keep each posed-frame inlier.
    Returns ``(pts, keep, res)`` over the full observation array; a one-line
    reprojection-inlier summary is printed as the only structure-quality signal
    (the pose-only save had none). A low inlier fraction means the poses/focal
    do not yield consistent multi-view structure even where the poses
    similarity-align to a reference."""
    pts = fill_new_points(pts, all_c, all_i, all_u, rvec, tvec, posed, f)
    ok = posed[all_i] & ~np.isnan(pts[:, 0])[all_c]
    res = np.full(len(all_c), np.inf)
    if ok.any():
        # (empty-input guard: ray_to_pixel_batch returns shape (0, 0) for an
        # empty batch, which does not broadcast against all_u[ok]'s (0, 2))
        xc = (
            Rotation.from_rotvec(rvec[all_i[ok]]).apply(pts[all_c[ok]])
            + tvec[all_i[ok]]
        )
        proj = make_cam(f).ray_to_pixel_batch(np.ascontiguousarray(xc))
        res[ok] = np.linalg.norm(proj - all_u[ok], axis=1)
    keep = res < TRIM_PX
    if not quiet:
        r = res[ok & np.isfinite(res)]
        if len(r):
            print(
                f"structure: {int(ok.sum())} dense obs triangulated; reproj "
                f"<2px {100 * (r < 2).mean():.1f}% <4px {100 * (r < 4).mean():.1f}% "
                f"<10px {100 * (r < 10).mean():.1f}% (median kept "
                f"{np.median(res[keep]):.2f} px)"
            )
        else:
            print("structure: 0 dense obs triangulated")
    return pts, keep, res


# ── Seed-stage snapshots (debug) ─────────────────────────────────────────────
# When SFMTOOL_SEED_SNAPSHOT_DIR is set, the SEED stage writes a numbered .sfmr
# at every checkpoint of its pipeline — stage 1's affine factorization, probe,
# widen and photometric verify (exp_fast_seed), the released estimate of each
# pass, and this module's photometric finalization (dense / embed / culled) — so
# every intermediate state can be opened in the SfM Explorer.  Unset (the
# default) every hook is a no-op and the run is byte-identical to production.
# Files are named `NN-<stage>[-<pass>-<attempt>].sfmr` with NN the checkpoint
# index, so a directory listing reads in pipeline order.


def seed_snapshot_path(tag):
    """Output path for a seed snapshot, or None when snapshots are disabled."""
    snap_dir = os.environ.get("SFMTOOL_SEED_SNAPSHOT_DIR")
    if not snap_dir:
        return None
    out = Path(snap_dir)
    out.mkdir(parents=True, exist_ok=True)
    return out / f"{tag}.sfmr"


def seed_snapshot(
    tag,
    data,
    f,
    rvec,
    tvec,
    pts,
    posed,
    extra_tool_options=None,
    path=None,
    release_grade=False,
):
    """Write one seed checkpoint as a dense, Explorer-loadable .sfmr.

    ``rvec``/``tvec``/``posed`` are indexed by ``data``'s image index and
    ``pts`` by its cluster index — a caller holding a thinned working set maps
    back to that frame first.  The structure is completed exactly as the final
    save does it (``dense_structure`` at the given poses, then ``save_sfmr``),
    so the file shows the state as the pipeline holds it.  Debug instrumentation
    must never kill the run it instruments, so every failure is caught and
    reported as a one-line warning.

    ``path`` names an explicit destination instead of the env-gated snapshot
    directory, for the artifacts a run writes unconditionally (the seed's
    per-hypothesis releases); ``release_grade`` writes poses and points only.
    """
    if path is None:
        path = seed_snapshot_path(tag)
    if path is None:
        return None
    try:
        global _CAM_WH
        if _CAM_WH is None:
            _CAM_WH = tuple(data["dims"][0])
        all_c, all_i, all_u = data["obs_c"], data["obs_i"], data["obs_uv"]
        pts_dense, keep, res = dense_structure(
            all_c,
            all_i,
            all_u,
            f,
            rvec,
            tvec,
            np.array(pts, dtype=np.float64, copy=True),
            posed,
            quiet=True,
        )
        # `save_sfmr` keeps points with >= 2 surviving observations; when none
        # does (typically a state with 0-1 posed frames) there is no
        # reconstruction to write, and the writer's empty arrays raise instead
        # of saying so.  The missing file IS the signal.
        n_alive = int((np.bincount(all_c[keep], minlength=len(pts_dense)) >= 2).sum())
        if n_alive == 0:
            print(
                f"  [seed-snapshot {tag}: {int(np.asarray(posed).sum())} posed, "
                f"{int(keep.sum())} obs, no multi-view points; skipped]"
            )
            return None
        save_sfmr(
            data,
            f,
            rvec,
            tvec,
            pts_dense,
            keep,
            res,
            path,
            tool_options=extra_tool_options,
            quiet=True,
            release_grade=release_grade,
        )
        print(
            f"  [seed-snapshot {path.name}: {int(np.asarray(posed).sum())} posed, "
            f"{len(np.unique(all_c[keep]))} pts, {int(keep.sum())} obs, f={f:.1f}]"
        )
        return path
    except Exception as exc:
        print(f"  [seed-snapshot {tag} FAILED: {type(exc).__name__}: {exc}]")
        return None


def triangulate(obs_c, obs_i, u, rot, trans, used, n_cl, f):
    """Ray-midpoint triangulation of every cluster from the posed images,
    via the batch triangulation binding (clusters with < 2 posed
    observations stay NaN)."""
    from sfmtool._sfmtool.analysis import triangulate_batch

    pts = np.full((n_cl, 3), np.nan)
    sel = used[obs_i]
    if not sel.any():
        return pts
    oc, oi, uv = obs_c[sel], obs_i[sel], u[sel]
    # World-space unit rays and camera centers: x_cam = R x + t, so the world
    # ray is Rᵀ·(canonical camera ray of the full pixel) and the center -Rᵀ t.
    d_loc = make_cam(f).pixel_to_ray_batch(np.ascontiguousarray(uv))
    dirs = np.einsum("nji,nj->ni", rot[oi], d_loc)
    centers = -np.einsum("nji,nj->ni", rot[oi], trans[oi])
    # obs_c is cluster-sorted, so the selection is CSR-ready.
    uniq, counts = np.unique(oc, return_counts=True)
    offsets = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
    result = triangulate_batch(
        np.ascontiguousarray(dirs), np.ascontiguousarray(centers), offsets
    )
    good = counts >= 2
    pts[uniq[good]] = np.asarray(result["points"])[good]
    return pts


# ── Save as .sfmr ────────────────────────────────────────────────────────────


def save_sfmr(
    data,
    f,
    rvec,
    tvec,
    pts,
    keep,
    res,
    out_path,
    return_alive=False,
    tool_options=None,
    quiet=False,
    release_grade=False,
    operation="cluster_bootstrap",
):
    """Write the bootstrap as an ``embedded_patches`` reconstruction.

    The bootstrap's observations are the cluster patches' *refined*
    positions, not the SIFT detections, so they are stored inline as
    ``keypoints_xy`` rather than as feature indexes into the ``.sift``
    files (which would silently resolve back to the unrefined seeds).

    ``tool_options`` merges extra entries into the file's metadata (a debug
    snapshot declares what state it holds); ``operation`` names the stage that
    produced the artifact, which a reader takes the file's provenance from, and
    stays the bootstrap's own unless the caller is a different stage;
    ``quiet`` suppresses the summary
    line so an instrumented run stays readable.  ``release_grade`` stops at
    poses and points — no keypoint passthrough and no patch frames — which is
    what an inspectable side artifact needs and is far cheaper: the surfel
    solve below opens every posed image's `.sift` affine array.
    """
    from sfmtool._sfmtool.reconstruction import SfmrReconstruction
    from sfmtool._workspace import load_workspace_config
    from sfmtool.colmap.convention import world_rotate_w
    from sfmtool.colmap.io import (
        _build_sfmr_data_dict,
        _resolve_workspace_and_sift,
        build_metadata,
        finite_positions_xyzw,
    )

    out_path = Path(out_path).resolve()
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
    keypoints_xy = data["obs_uv"][ko].astype(np.float32)
    point_idx = np.array([remap[int(c)] for c in obs_c[ko]])
    obs_counts = np.bincount(point_idx, minlength=len(alive))

    # Write only the posed cameras. Every kept observation belongs to a posed
    # frame, so the images that carry one are exactly the posed set; the rest
    # hold the shared default seed pose (identical camera centers), which both
    # misrepresents the reconstruction (phantom cameras at the origin) and, when
    # many frames are unposed (e.g. an early growth snapshot), overflows the
    # viewer's k-d tree spatial index. Compact the per-image arrays and remap
    # the observation image indexes; everything downstream then works in the
    # posed-only image space (embedded_patches cannot be image-subset later).
    posed_img = np.unique(track_img)
    if len(posed_img) < len(names):
        img_remap = np.full(len(names), -1, np.int64)
        img_remap[posed_img] = np.arange(len(posed_img))
        names = [names[j] for j in posed_img]
        rvec = rvec[posed_img]
        tvec = tvec[posed_img]
        track_img = img_remap[track_img].astype(track_img.dtype)

    positions = pts[alive]
    per_point_err = np.zeros(len(alive), dtype=np.float32)
    np.add.at(per_point_err, point_idx, res[ko].astype(np.float32))
    per_point_err /= np.maximum(obs_counts, 1)

    # The internal poses are already canonical camera frame, in the COLMAP-world
    # gauge; only the world rotation W remains to reach the .sfmr canonical
    # world.  W rotates the point positions and, applied to each rotation row,
    # right-multiplies the world->camera rotations (R_int·Wᵀ); the camera-frame
    # translation is unchanged.
    rot_int = Rotation.from_rotvec(rvec).as_matrix()
    q_can = Rotation.from_matrix(
        world_rotate_w(rot_int.reshape(-1, 3)).reshape(-1, 3, 3)
    ).as_quat()[:, [3, 0, 1, 2]]
    t_can = tvec
    p_can = world_rotate_w(positions)

    (
        workspace_dir,
        _contents,
        resolved_names,
        ft_hashes,
        sc_hashes,
        thumbnails,
    ) = _resolve_workspace_and_sift(names, WS.resolve())

    # Colors from the .sift thumbnails at the (scaled) observation position.
    colors = np.zeros((len(alive), 3), dtype=np.uint8)
    uv = data["obs_uv"][ko]
    for k in range(len(ko)):
        th = np.asarray(thumbnails[track_img[k]])
        ty = int(np.clip(uv[k, 1] * th.shape[0] / h, 0, th.shape[0] - 1))
        tx = int(np.clip(uv[k, 0] * th.shape[1] / w, 0, th.shape[1] - 1))
        colors[point_idx[k]] = th[ty, tx]

    # The context camera, NOT a hardcoded SIMPLE_PINHOLE.  The structure this
    # writes was triangulated through ``make_cam`` (equidistant rays under a
    # fisheye context), so stamping a pinhole model on it hands every
    # downstream consumer — the photometric embed, the reprojection culls, the
    # finalization BA — a camera that does not describe the observations.
    # Identical to the previous literal on the pinhole default.
    camera = make_cam(float(f))

    opts = {"camera_model": _CAM_CONTEXT["model"], "focal_grid": F_GRID}
    opts.update(tool_options or {})
    metadata = build_metadata(
        workspace_dir=workspace_dir,
        output_path=out_path,
        workspace_config=load_workspace_config(workspace_dir),
        operation=operation,
        tool_name="sfmtool",
        tool_options=opts,
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
        # A release-grade file is read on its own, often against a workspace
        # whose `.sift` files were never written (the observations came from the
        # cluster patches). Carry the refined positions inline so its
        # observations stay readable; the patch-bearing path below states them
        # as an `embedded_patches` file instead.
        keypoints_xy=keypoints_xy if release_grade else None,
        point_indexes=point_idx,
        observation_counts=obs_counts,
        feature_tool_hashes=ft_hashes,
        sift_content_hashes=sc_hashes,
        thumbnails=thumbnails,
        metadata=metadata,
    )

    recon = SfmrReconstruction.from_data(workspace_dir, sfmr_dict)

    if release_grade:
        # Poses and points, and stop: no keypoint passthrough, no patch cloud,
        # no normals.  The surfel solve below is the expensive half of this
        # writer (a `.sift` read per posed image plus a per-point least
        # squares), and a side artifact meant for inspection carries none of it.
        recon.save(out_path)
        if not quiet:
            print(f"\nwrote {out_path} ({len(alive)} points, release-grade)")
        return (recon, alive, posed_img) if return_alive else recon

    # ── Surfel frames copied from the cluster patches ────────────────────
    # Each member's stored 2x2 warp is the projection of the cluster's
    # common surfel into that image, so the 3D patch frame is recoverable:
    # solve J_k·B = A_k per point (J_k the projection Jacobian at the
    # point, B the 3x2 map from reference-image pixels to 3D on the surfel
    # plane; the reference row contributes J_ref·B = I), then
    # u = B·(r, 0), v = B·(0, r) with r the refinement radius in reference
    # pixels (keypoint-frame radius x the reference feature's scale).
    from sfmtool._sfmtool.patches import PatchCloud
    from sfmtool._sfmtool.io import read_sift, read_sift_metadata
    from sfmtool.colmap.convention import world_rotate_w
    from sfmtool.sift.file import get_sift_path_for_image

    feature_scales = {}
    image_file_hashes = []
    for i, name in enumerate(names):
        sp = get_sift_path_for_image(workspace_dir / name)
        meta = read_sift_metadata(sp)["metadata"]
        image_file_hashes.append(bytes.fromhex(meta["image_file_xxh128"]))
        shapes = np.asarray(read_sift(sp)["affine_shapes"], dtype=np.float64)
        feature_scales[i] = 0.5 * (
            np.linalg.norm(shapes[:, :, 0], axis=1)
            + np.linalg.norm(shapes[:, :, 1], axis=1)
        )

    # The surfel-frame solve below runs in the COLMAP +Z-forward camera frame,
    # so it runs on the poses flipped back to that frame by S = diag(1, -1, -1);
    # its world-space u/v/normal outputs convert to the canonical world by the
    # same W as the points, at the end.  Positions stay in the COLMAP-world
    # gauge.
    #
    # There are two arms.  The PERSPECTIVE arm is the historical one, written
    # against the pinhole projection Jacobian `(f/z)[I | -p_r]` and its
    # `z > 0` in-front test, with the reference right-inverse `(z/f)[I; 0]` and
    # the null direction `X/z`.  Under a RAY-PATH context (fisheye) both are
    # wrong: `z` is not a distance past 90 deg off axis (it crosses zero and
    # then goes negative over the very periphery a >180 deg capture exists to
    # image), and the pinhole Jacobian describes a different map.  The fisheye
    # arm is the same solve stated model-generically -- the camera's own 2x3
    # Jacobian, its minimum-norm right inverse, and the unit viewing ray as the
    # null direction -- which reduces to the perspective arm up to a
    # reparameterization of the out-of-plane slope.  Both are kept because the
    # Tikhonov prior and the obliquity cap live on that slope, so the pinhole
    # path stays bit-identical only if its own parameterization does.
    fisheye_frames = fisheye_stage1()
    s_flip = np.array([1.0, -1.0, -1.0])
    rot_all = Rotation.from_rotvec(rvec).as_matrix() * s_flip[None, :, None]
    tvec_col = tvec * s_flip
    warps = data["obs_warp"][ko]
    is_ref = data["obs_ref"][ko]
    radius_kf = data["refine_radius"]
    half_u = np.zeros((len(alive), 3), dtype=np.float64)
    half_v = np.zeros((len(alive), 3), dtype=np.float64)
    normals = np.zeros((len(alive), 3), dtype=np.float64)
    p_starts = np.searchsorted(point_idx, np.arange(len(alive) + 1))
    # The reference constraint J_ref·B = I determines B up to a 2-vector
    # b_z — the surfel's out-of-plane slope in the reference camera frame
    # (B = R_refᵀ·[(z_r/f)·I + p_r·b_z ; b_z] with p_r the normalized ref
    # coords).  Each other member contributes A_k − (z_r/f)·M2 = c_k·b_z
    # with M = J_k·R_refᵀ = [M2 | m3] and c_k = M2·p_r + m3.  The tilt is
    # exactly the depth-like weakly-observed direction, so the solve gets a
    # fronto-parallel Tikhonov prior (weight relative to the members'
    # leverage) and a hard obliquity cap; these are what the photometric
    # normal refinement later polishes.
    tan_cap = np.tan(np.radians(80.0))
    # Under the fisheye context, every kept observation's camera-frame point and
    # the camera model's own 2x3 pixel Jacobian there, taken ONCE for the whole
    # writer by central difference of ``ray_to_pixel`` — the same measure
    # ``WarpMap`` takes its warp Jacobian by, so the frames agree with the render.
    xc_obs = j_obs = None
    if fisheye_frames:
        xc_obs = (
            np.einsum("nij,nj->ni", rot_all[track_img], positions[point_idx])
            + tvec_col[track_img]
        )
        j_obs = _colmap_proj_jacobian(make_cam(float(f)), xc_obs, s_flip)
    # Patch EXTENT for every point at once, through THE CAMERA'S OWN sizing rule
    # (``CameraIntrinsics.pixel_radius_to_world_batch``): the world size that
    # subtends the detection's pixel radius at the reference view.  One rule for
    # every model — it reduces to `r_px*|z|/f` for a pinhole and `r_px*d/f` under
    # the equidistant map — and, more to the point, it is the SAME
    # implementation ``PatchExtent::FeatureSize`` uses, so this writer and the
    # core cannot drift apart.  (They did once: the writer carried a hand-copied
    # pinhole formula through Phase 3b, which is exactly how a fisheye seed came
    # to be sized through a pinhole camera.)
    ref_row = np.full(len(alive), -1, np.int64)
    for p in range(len(alive)):
        lo, hi = int(p_starts[p]), int(p_starts[p + 1])
        here = np.nonzero(is_ref[lo:hi])[0]
        if len(here):
            ref_row[p] = lo + int(here[0])
    ext_all = np.zeros(len(alive))
    has_ref = ref_row >= 0
    if has_ref.any():
        rr = ref_row[has_ref]
        ri = track_img[rr]
        xc_ref_all = (
            np.einsum("nij,nj->ni", rot_all[ri], positions[has_ref]) + tvec_col[ri]
        )
        r_px_all = np.array(
            [
                radius_kf * feature_scales[int(i)][int(track_feat[k])]
                for i, k in zip(ri.tolist(), rr.tolist())
            ]
        )
        ext_all[has_ref] = make_cam(float(f)).pixel_radius_to_world_batch(
            # The camera projects CANONICAL rays; the writer works in the COLMAP
            # frame, and the two differ by the involution S = diag(1, -1, -1).
            np.ascontiguousarray(xc_ref_all * s_flip),
            np.ascontiguousarray(r_px_all),
        )
    for p in range(len(alive)):
        lo, hi = int(p_starts[p]), int(p_starts[p + 1])
        refs_here = np.nonzero(is_ref[lo:hi])[0]
        if len(refs_here) == 0:
            continue  # reference member trimmed: leave the zero (no-patch) frame
        k_ref = lo + int(refs_here[0])
        i_ref = int(track_img[k_ref])
        x_pt = positions[p]
        r_ref = rot_all[i_ref]
        xc_ref = r_ref @ x_pt + tvec_col[i_ref]
        if fisheye_frames:
            # RANGE, not optical-axis depth: the distance an angular size
            # ``r_px / f`` actually spans, and the one quantity that stays
            # positive over a >180 deg field.
            d_ref = max(float(np.linalg.norm(xc_ref)), 1e-6)
            n_ref = xc_ref / d_ref  # unit viewing ray == null direction of J_ref
            scale_ref = d_ref / f
            b_perp = np.linalg.pinv(j_obs[k_ref])  # 3x2, columns ⊥ n_ref
        else:
            z_ref = max(xc_ref[2], 1e-6)
            p_r = xc_ref[:2] / z_ref
            scale_ref = z_ref / f
        rows, rhs = [], []
        for k in range(lo, hi):
            if k == k_ref:
                continue
            i = int(track_img[k])
            if fisheye_frames:
                m = j_obs[k] @ rot_all[i] @ r_ref.T
                c_k = m @ n_ref
                resid = warps[k] - m @ b_perp
            else:
                xc = rot_all[i] @ x_pt + tvec_col[i]
                z = max(xc[2], 1e-6)
                j_proj = (f / z) * np.array(
                    [[1.0, 0.0, -xc[0] / z], [0.0, 1.0, -xc[1] / z]]
                )
                m = j_proj @ rot_all[i] @ r_ref.T
                c_k = m[:, :2] @ p_r + m[:, 2]
                resid = warps[k] - scale_ref * m[:, :2]
            for j in range(2):
                rows.append([c_k[0] * (1 - j), c_k[0] * j])
                rows.append([c_k[1] * (1 - j), c_k[1] * j])
                rhs.append(resid[0, j])
                rhs.append(resid[1, j])
        if not rows:
            continue
        rows = np.asarray(rows)
        rhs = np.asarray(rhs)
        # Fronto prior: damping rows scaled to a fraction of member leverage.
        lam = WRITER_FRONTO_LAM * np.sqrt((rows**2).sum() / max(len(rows), 1))
        rows = np.vstack([rows, [[lam, 0.0], [0.0, lam]]])
        rhs = np.concatenate([rhs, [0.0, 0.0]])
        b_z = np.linalg.lstsq(rows, rhs, rcond=None)[0]
        # Obliquity cap: tan(tilt) = |b_z| / (the in-plane scale of the
        # reference right-inverse) — z_r/f for the pinhole, d_r/f (the radial
        # singular value of the equidistant map's inverse) for the fisheye arm.
        b_norm = np.linalg.norm(b_z)
        max_bz = tan_cap * scale_ref
        if b_norm > max_bz:
            b_z *= max_bz / b_norm
        if fisheye_frames:
            b_map = r_ref.T @ (b_perp + np.outer(n_ref, b_z))
        else:
            b_map = r_ref.T @ np.vstack(
                [scale_ref * np.eye(2) + np.outer(p_r, b_z), b_z[None, :]]
            )
        # The tilt solve above determines the surfel's NORMAL and nothing else.
        # The frame itself is a SQUARE: a SIFT detection is a round region of the
        # image, so the surface element it seeds is a square patch that
        # orientation may TILT but nothing may DISTORT.  Foreshortening is the
        # projection's business — it belongs in the render, not in the stored
        # frame.  (b_map's own two columns are a sheared, anisotropic
        # parallelogram: they carry the tilt's foreshortening inside the extents.
        # Storing those made the seed's patches stretched — median |u|/|v| 2.6 on
        # 20240614_224422531 with 22% beyond 4x — and worse, `refine_normals`
        # later re-solves the normal to near-fronto while KEEPING the extents, so
        # the patch ends up facing the camera and still stretched: the measured
        # tilt is 0.8-7.6 deg median across the fleet while the elongation stays.)
        n3 = np.cross(b_map[:, 0], b_map[:, 1])
        norm = np.linalg.norm(n3)
        if norm < 1e-12:
            continue
        n3 /= norm
        cam_c = -r_ref.T @ tvec_col[i_ref]
        if np.dot(n3, cam_c - x_pt) < 0:
            n3 = -n3  # front-facing; the frame below is built around the normal
        # In-plane u: the reference image's x direction projected onto the plane.
        # This keeps the upright-bitmap convention structural rather than
        # corrective — bitmap columns run along reference-image x, and
        # v = n x u then comes out along -image-y, which is exactly what the
        # raster's row reversal (rows step along -v) expects.
        ax = r_ref.T @ np.array([1.0, 0.0, 0.0])
        u3 = ax - np.dot(ax, n3) * n3
        nu = np.linalg.norm(u3)
        if nu < 1e-9:  # plane edge-on to image x: fall back to image y
            ay = r_ref.T @ np.array([0.0, 1.0, 0.0])
            u3 = ay - np.dot(ay, n3) * n3
            nu = np.linalg.norm(u3)
            if nu < 1e-9:
                continue
        u3 /= nu
        v3 = np.cross(n3, u3)  # u x v == n exactly, for a unit u perpendicular to n
        # ONE scalar extent: the detection's fronto world size at the reference
        # view, taken from the camera itself in the batch above (`ext_all`).
        # b_map's own fronto case (b_z = 0) is r_ref^T (z/f) I, so a detection of
        # r_px reference pixels spans exactly the camera's own
        # `pixel_radius_to_world` at that point — the same detection-time
        # quantity the duplicate collapse sizes its radius from, now the surfel's
        # only size parameter.  Under the equidistant map that is the RANGE form
        # (r_px pixels subtend r_px / f radians everywhere, spanning
        # r_px * d_ref / f at range d_ref); the optical-axis form is that times
        # cos(theta_ref), which shrinks the patch two-fold at 60 deg off axis and
        # to nothing at 90 — a zero-extent frame is written as NO patch at all
        # (`PatchCloud.from_halfvec_arrays` drops a zero `u` row).
        ext = float(ext_all[p])
        half_u[p], half_v[p], normals[p] = u3 * ext, v3 * ext, n3

    # COLMAP -> canonical for the direction quantities (same W as the points).
    half_u = np.asarray(world_rotate_w(half_u), dtype=np.float32)
    half_v = np.asarray(world_rotate_w(half_v), dtype=np.float32)
    normals = np.asarray(world_rotate_w(normals), dtype=np.float32)

    cloud = PatchCloud.from_halfvec_arrays(half_u, half_v, np.asarray(p_can))
    recon = recon.clone_with_changes(
        feature_source="embedded_patches",
        keypoints_xy=keypoints_xy,
        image_file_hashes=image_file_hashes,
        normals=normals,
        patches=cloud,
    )
    recon.save(out_path)
    n_patched = int(np.count_nonzero(np.linalg.norm(half_u, axis=1) > 0))
    if not quiet:
        print(
            f"\nwrote {out_path} ({len(alive)} points, {int(obs_counts.sum())} obs, "
            f"{recon.feature_source}, {n_patched} warp-derived patch frames)"
        )
    # `alive` maps output point index -> source cluster index; `posed_img` maps
    # output image index -> data image index (the writer keeps only posed
    # frames). The seed finalization maps refined points/views back through them.
    return (recon, alive, posed_img) if return_alive else recon
