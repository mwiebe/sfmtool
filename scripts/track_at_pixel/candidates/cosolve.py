# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Cosolve: ``renormal``'s track, re-oriented over a neighbourhood of fitted patches.

The prototypes of ``specs/drafts/surface-co-solve.md`` that need more than one
patch, applied in this order to the track ``renormal`` returns
(``renormal_options`` are passed to it, and ``renormal_source`` sets its
normal sources). The defaults are the combination the harness kept:
``renormal`` with its neighbours admitted by position (prototype 1), and the
seeded grid when the reconstruction has no point near the track.

- ``size="curve"`` (prototype 5, not kept): the track is resized to each of
  ``size_steps`` times its half-size, given an anchored fit and read; the
  largest size whose median ZNCC is within ``size_drop`` of the best and that
  passes the gates is kept. While the patch stays on one surface the ZNCC
  holds as it grows; once it crosses an edge, the part beyond moves
  differently between views and the ZNCC drops.
- ``split_source`` (prototype 4 alone, not kept): the normal from split
  patches along both in-plane axes. Two copies of half the size, half a size
  either side of the centre along an axis, are fitted; the line between
  their fitted centres lies in the surface, and two such lines fix the
  normal. ``split_rounds`` rounds of estimate, tilt and anchored fit.
- ``grid`` (prototype 3, the seeded grid): ``"3x3"`` or ``"cross"``, run
  according to ``grid_when``. Copies of the track shifted ``grid_spacing``
  half-sizes (half a diameter) along its in-plane axes are fitted. Each round
  the normal comes from all the fitted centres (a plane through them; where
  they lie along a line, only the turn that line determines), and the
  patches are tilted to it and fitted again, the centre anchored on the
  pixel. A copy whose fit fails or lands more than ``grid_outlier_halves``
  half-sizes off the plane of the others is dropped: the surface ends there.
  With ``grid_split`` (prototype 4 in the grid, not kept), an axis the
  centres leave free is fixed by a split patch across it. With
  ``grid_recon_points``, the reconstruction's own points near the track join
  the grid by their positions only, when they lie within
  ``grid_recon_plane_halves`` of the track's plane. After the re-orientation
  (``regrow``) the geometry search runs again from the query: a pose close to
  right finds views a wrong one missed.

The result replaces ``renormal``'s track only when it passes the same gates
(and, with ``min_zncc_gain``, reads well enough against it).
"""

from __future__ import annotations

import numpy as np

from api import TrackAtPixelResult
from candidates import renormal
from candidates.common import (
    anchored_fit,
    clean,
    median_zncc,
    patch_normal,
    shape_track,
)

DEFAULTS = {
    "renormal_options": {},
    # Shorthand for renormal_options["source"].
    "renormal_source": "nbpos+nb3dpos+photo",
    "size": "none",
    "size_steps": [1.0, 1.5, 2.0, 3.0],
    "size_drop": 0.02,
    "split_source": False,
    "split_rounds": 2,
    "grid": "3x3",
    # "always"; "photo": only when renormal's normal came from photometry (no
    # reconstructed neighbour gave one); "empty": only when the reconstruction
    # holds no point at all to lean on, as early in a reconstruction.
    "grid_when": "empty",
    "grid_spacing": 1.0,
    "grid_rounds": 2,
    "grid_recon_radius": 2.0,
    "grid_outlier_halves": 0.75,
    "grid_split": False,
    # Refuse the grid's normal when its fitted centres lie further than this
    # (rms, in spacings) off their plane; None keeps it always.
    "grid_max_resid": None,
    "grid_recon_points": False,
    "grid_recon_plane_halves": 0.5,
    # Spread a direction of the fitted centres needs, in half-sizes, before it
    # counts as fixing the normal along it.
    "fixed_spread_halves": 0.25,
    "stop_deg": 1.0,
    # After the re-orientation, the geometry search from the query and an
    # anchored refit: a pose close to right finds views a wrong one missed.
    "regrow": True,
    # When set, the result replaces renormal's track only when its median ZNCC
    # is at least renormal's plus this (negative: a tolerance).
    "min_zncc_gain": None,
}


def _unit(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else None


def _face(n, to_cam):
    return n if n @ to_cam >= 0 else -n


def _half(track) -> float:
    return float(np.linalg.norm(track.placement["u_halfvec"]))


def _axes(track):
    pl = track.placement
    return _unit(pl["u_halfvec"]), _unit(pl["v_halfvec"])


def _move(ctx, track, d):
    from sfmtool._sfmtool import bench as B

    u, v = _axes(track)
    m = _unit(np.cross(u, v))
    track, _ = B.translate_patch(
        track, ctx.edited, (float(d @ u), float(d @ v), float(d @ m))
    )
    return track


def _fit(ctx, track):
    """A free fit, or ``None`` when it is refused or leaves under two views."""
    from sfmtool._sfmtool import bench as B

    try:
        track, _ = B.fit(track, ctx.edited, ctx.pyramids)
        track, _ = B.apply_thresholds(track)
    except ValueError:
        return None
    if track.at_infinity or track.verdict_counts[0] < 2:
        return None
    return track


def _tilt(ctx, track, n):
    from sfmtool._sfmtool import bench as B

    try:
        track, _ = B.tilt_patch(track, ctx.edited, tuple(float(x) for x in n))
    except ValueError:
        pass
    return track


def _deg(a, b):
    return float(np.degrees(np.arccos(np.clip(abs(float(a @ b)), -1, 1))))


def plane_normal(points, current, half, opts):
    """The normal a set of fitted centres fixes, keeping what they leave free.

    Returns ``(normal, fixed)``: the in-plane directions whose spread is at
    least ``fixed_spread_halves`` half-sizes. Two fix the normal; one fixes
    only the turn about its perpendicular, and the tilt about the line itself
    is left as ``current`` had it.
    """
    P = np.asarray(points) - np.mean(points, axis=0)
    vals, vecs = np.linalg.eigh(P.T @ P / len(points))
    fixed = [
        vecs[:, k]
        for k in range(3)
        if np.sqrt(max(vals[k], 0.0)) > opts["fixed_spread_halves"] * half
    ]
    if len(fixed) >= 2:
        n = _unit(np.cross(fixed[-1], fixed[-2]))
    elif len(fixed) == 1:
        a = fixed[0]
        n = _unit(current - (current @ a) * a)
    else:
        return current, []
    if n is None:
        return current, []
    return (n if n @ current > 0 else -n), fixed[-2:]


def split_line(ctx, track, axis):
    """The in-surface direction across ``axis``, from two fitted half-size copies."""
    from sfmtool._sfmtool import bench as B

    half = _half(track)
    centres = []
    for sign in (1.0, -1.0):
        try:
            s, _ = B.resize_patch(track, ctx.edited, half / 2)
        except ValueError:
            return None
        s = _fit(ctx, _move(ctx, s, sign * (half / 2) * axis))
        if s is None:
            return None
        centres.append(np.asarray(s.position, float))
    return _unit(centres[0] - centres[1])


def split_normal(ctx, track, n):
    """The normal from split patches along both in-plane axes."""
    u, v = _axes(track)
    lines = [split_line(ctx, track, a) for a in (u, v)]
    lines = [b for b in lines if b is not None]
    if len(lines) == 2:
        m = _unit(np.cross(lines[0], lines[1]))
    elif len(lines) == 1:
        b = lines[0]
        m = _unit(n - (n @ b) * b)
    else:
        return n
    if m is None:
        return n
    return m if m @ n > 0 else -m


def gated(track, q, pixel, ropts) -> bool:
    return renormal.passes_gates(track, q, pixel, ropts)


def size_curve(ctx, track, q, pixel, opts, ropts, info):
    from sfmtool._sfmtool import bench as B

    h0 = _half(track)
    tried = []
    for s in opts["size_steps"]:
        try:
            t = track if s == 1.0 else B.resize_patch(track, ctx.edited, h0 * s)[0]
            t = anchored_fit(ctx, t, q, pixel, ropts["anchor_refits"])
            t, _ = B.apply_thresholds(t)
        except ValueError:
            continue
        tried.append((s, median_zncc(t), gated(t, q, pixel, ropts), t))
    ok = [x for x in tried if x[2]]
    info["size_curve"] = [[s, z, g] for s, z, g, _ in tried]
    if not ok:
        return track
    zbest = max(z for _, z, _, _ in ok)
    s, z, _, t = max(
        (x for x in ok if x[1] >= zbest - opts["size_drop"]), key=lambda x: x[0]
    )
    info["size_scale"] = s
    return t


def grid_normal(ctx, track, q, pixel, image, opts, ropts, info):
    """The seeded grid: the centre's normal after the repeated step over fitted copies."""
    cam = ctx.camera(image)
    half = _half(track)
    spacing = opts["grid_spacing"] * half
    u, v = _axes(track)
    offsets = (
        [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1) if (i, j) != (0, 0)]
        if opts["grid"] == "3x3"
        else [(1, 0), (-1, 0), (0, 1), (0, -1)]
    )
    seeds = [
        _fit(ctx, _move(ctx, track, spacing * (i * u + j * v))) for i, j in offsets
    ]
    patches = [track] + [s for s in seeds if s is not None]
    info["grid_seeded"] = len(patches) - 1
    extra = []
    if opts["grid_recon_points"]:
        n0 = patch_normal(track)
        X0 = np.asarray(track.position, float)
        for o in ctx.points_near(X0, k=12, radius=opts["grid_recon_radius"] * spacing):
            if (
                abs(float((o["position"] - X0) @ n0))
                <= opts["grid_recon_plane_halves"] * half
            ):
                extra.append(np.asarray(o["position"], float))
        info["grid_recon_points"] = len(extra)
    centre_normals = []
    for rnd in range(opts["grid_rounds"]):
        X = np.asarray([p.position for p in patches], float)
        to_cam = _unit(cam.center - X[0])
        normals = [_face(patch_normal(p), to_cam) for p in patches]
        # Drop copies far off the plane of the others: the surface ends there.
        if len(patches) >= 5:
            n_all, _ = plane_normal(X, normals[0], half, opts)
            resid = np.abs((X - X.mean(0)) @ n_all)
            keep = [0] + [
                k
                for k in range(1, len(patches))
                if resid[k] <= opts["grid_outlier_halves"] * half
            ]
            if len(keep) < len(patches):
                info.setdefault("grid_dropped", 0)
                info["grid_dropped"] += len(patches) - len(keep)
                patches = [patches[k] for k in keep]
                X, normals = X[keep], [normals[k] for k in keep]
        # The grid is small enough to be one plane: every patch takes the
        # normal of all the fitted centres. A copy's fit may slide far along
        # the surface (a two-view track on the ground slides by several
        # half-sizes) and still lie in the plane, which is all that is read.
        pts = list(X) + extra
        n, fixed = plane_normal(pts, normals[0], half, opts)
        free_axis = _unit(np.cross(n, fixed[0])) if len(fixed) == 1 else None
        if len(pts) < 3:
            n = normals[0]
        new = [n if n @ m > 0 else -n for m in normals]
        if opts["grid_split"] and free_axis is not None:
            b = split_line(ctx, patches[0], free_axis)
            if b is not None:
                a = _unit(np.cross(new[0], free_axis))
                m = _unit(np.cross(a, b))
                if m is not None:
                    new = [m if m @ k > 0 else -m for k in new]
                    info["grid_split_used"] = info.get("grid_split_used", 0) + 1
        turned = max(_deg(a, b) for a, b in zip(new, normals))
        centre_normals.append(new[0])
        refit = []
        for k, (p, n) in enumerate(zip(patches, new)):
            t = _tilt(ctx, p, n)
            if k == 0:
                try:
                    t = anchored_fit(ctx, t, q, pixel, ropts["anchor_refits"])
                except ValueError:
                    t = p
            else:
                t = _fit(ctx, t)
            if t is not None:
                refit.append(t)
            elif k == 0:
                refit.append(p)
        patches = refit
        if turned < opts["stop_deg"]:
            break
    info["grid_rounds_run"] = len(centre_normals)
    info["grid_kept"] = len(patches) - 1
    if not centre_normals:
        return None
    # How far the final fitted centres lie off their plane, in spacings. On a
    # surface they lie on it; on something with no surface (a distant chimney)
    # they scatter along the rays and the plane through them means nothing.
    X = np.asarray([p.position for p in patches], float)
    n = centre_normals[-1]
    resid = float(np.sqrt(np.mean(((X - X.mean(0)) @ n) ** 2))) / spacing
    info["grid_resid"] = resid
    if opts["grid_max_resid"] is not None and (
        len(patches) < 4 or resid > opts["grid_max_resid"]
    ):
        info["grid_refused"] = True
        return None
    return n


def build_track(ctx, image: int, pixel, options: dict | None = None):
    from sfmtool._sfmtool import bench as B

    opts = {**DEFAULTS, **(options or {})}
    given = dict(opts["renormal_options"])
    if opts["renormal_source"]:
        given["source"] = opts["renormal_source"]
    ropts = {**renormal.DEFAULTS, **given}
    result = renormal.build_track(ctx, image, pixel, given)
    track, q = result.track, result.query_observation
    if track.at_infinity or track.placement is None:
        return result
    info: dict = {}
    result.diagnostics["cosolve"] = info
    start = track
    try:
        if opts["size"] == "curve":
            track = size_curve(ctx, track, q, pixel, opts, ropts, info)
        normal = None
        if opts["split_source"]:
            n = patch_normal(track)
            for _ in range(opts["split_rounds"]):
                n = split_normal(ctx, track, n)
                track = shape_track(ctx, track, q, pixel, n)
                track = anchored_fit(ctx, track, q, pixel, ropts["anchor_refits"])
            normal = n
        source = (result.diagnostics.get("renormal") or {}).get("source")
        info["renormal_source"] = source
        when = opts["grid_when"]
        if opts["grid"] != "none" and (
            when == "always"
            or (when == "photo" and source in ("photo", "photowide", None))
            or (when == "empty" and not ctx.points_near(track.position, k=1))
        ):
            normal = grid_normal(ctx, track, q, pixel, image, opts, ropts, info)
        if normal is not None:
            before = patch_normal(track)
            info["turn_deg"] = _deg(normal, before)
            track = shape_track(ctx, track, q, pixel, normal)
            track = anchored_fit(ctx, track, q, pixel, ropts["anchor_refits"])
            track, _ = B.apply_thresholds(track)
            if opts["regrow"] and track.verdict_counts[0] >= 2:
                grown, geo = B.search_geometry(track, q, ctx.edited, ctx.pyramids)
                info["regrow_added"] = geo["added"]
                if geo["added"]:
                    grown, _ = B.evaluate(grown, ctx.edited, ctx.pyramids)
                    grown, _ = B.apply_thresholds(grown)
                    grown = anchored_fit(ctx, grown, q, pixel, ropts["anchor_refits"])
                    track, _ = B.apply_thresholds(grown)
        if track is not start and ropts["clean_after"]:
            track = clean(ctx, track, q, pixel, ropts, info)
    except ValueError as e:
        info["error"] = str(e)
        return result
    ok = gated(track, q, pixel, ropts)
    info["gates"] = ok
    if not ok:
        return result
    z0, z1 = median_zncc(start), median_zncc(track)
    info.update(zncc_before=z0, zncc_after=z1)
    if opts["min_zncc_gain"] is not None and z1 < z0 + opts["min_zncc_gain"]:
        info["kept"] = False
        return result
    return TrackAtPixelResult(
        track=track, query_observation=q, diagnostics=result.diagnostics
    )


__all__ = ["build_track", "DEFAULTS"]
