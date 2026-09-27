# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Staged: the position first, then the views, then the normal.

Prototype 9 of ``specs/drafts/surface-co-solve.md``, the order that worked on
the query traced in its "One query, step by step". The cascade's track (with
``renormal``'s fallbacks) is taken for its position only, and then:

1. **Vetted clusters** (``vet_clusters``). The clusters with a member in the
   queried image within ``vet_radius_px`` of the pixel are checked with the
   posed cameras: the queried image's member must be the reference or kept,
   and the reference and kept members must triangulate in front of every
   camera with every reprojection error under ``vet_max_reproj_px`` (the worst
   member is dropped and the rest tried again while three or more remain). The
   pixel is carried through each cluster that passes into the images the
   track does not see yet, as candidates the reading then judges.
2. **Growth** (``grow_rounds``). The geometry search from the queried
   sighting, a reading and the thresholds, then the **relaxed** views turned
   in, then an anchored fit; repeated while the search adds views. A view is
   relaxed in when the thresholds turned it out but it reads at least
   ``relax_min_zncc``, its correlation peak sits within
   ``relax_max_shift_px`` of its keypoint and its keypoint within
   ``relax_max_projection_px`` of the point's projection: it lands where the
   geometry says, and a ZNCC under the bar is what a wrong normal does to a
   view that sees the patch from a different direction.
3. **The normal** (``normal``, a chain as in ``renormal``). ``nbpos`` and
   ``nb3dpos`` take the reconstruction's neighbours by position (prototype 1).
   ``photo_grid`` is the photometric search on the grown track, then the
   seeded grid laid out on its estimate (``grid_rounds`` rounds; with
   ``grid_when="empty"`` only when no reconstructed point is near the
   track). The photometric and grid steps are kept only when the track's
   reading does not fall: the median ZNCC less ``accept_shift_weight`` times
   the median correlation peak offset, with ``accept_tolerance``. The peak
   offset separates normals the ZNCC cannot. A neighbours' normal is kept
   whatever the reading says (``accept_neighbours="always"``): on a track with
   few views the right normal can read lower than a wrong one.
4. **Growth again** (``regrow_rounds``), now that the patch faces the right
   way, then the relaxed views judged again at ``final_min_zncc``, the
   cleaning and ``renormal``'s gates.

A track that fails the gates, or a step that raises, falls back to
``cosolve``'s default result.
"""

from __future__ import annotations

import numpy as np

from api import TrackAtPixelResult
from candidates import cosolve, renormal
from candidates.common import (
    anchored_fit,
    clean,
    in_frame,
    median_zncc,
    shape_track,
)

DEFAULTS = {
    "vet_clusters": True,
    "vet_radius_px": 24.0,
    "vet_max_clusters": 8,
    "vet_max_reproj_px": 2.0,
    "grow_rounds": 2,
    "relax": True,
    "relax_min_zncc": 0.72,
    "relax_max_shift_px": 1.5,
    "relax_max_projection_px": 1.5,
    "normal": "nbpos+nb3dpos+photo_grid",
    "grid_rounds": 2,
    # "always", or "empty": only when the reconstruction holds no point near
    # the track (in the full pass the grid loses to single-patch photometry).
    "grid_when": "empty",
    # "always", or "reading": a neighbours' normal only when the reading
    # accepts it, as the photometric and grid normals are.
    "accept_neighbours": "always",
    "accept_shift_weight": 0.02,
    "accept_tolerance": 0.005,
    "regrow_rounds": 1,
    "final_min_zncc": 0.8,
}


def _unit(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else None


def _row(o, key):
    v = (o.get("track") or {}).get(key)
    return None if v is None else float(v)


def reading(track, weight) -> float:
    """Median ZNCC over the ``in`` views less ``weight`` times their median peak offset."""
    z = median_zncc(track)
    shifts = [
        _row(o, "seed_shift_px")
        for o in track.observations
        if o["verdict"] == "in" and _row(o, "seed_shift_px") is not None
    ]
    return z - weight * (float(np.median(shifts)) if shifts else 0.0)


def triangulate(ctx, sightings):
    """The least-squares meeting point of ``(image, pixel)`` rays, and each reprojection error."""
    A, b = np.zeros((3, 3)), np.zeros(3)
    for image, px in sightings:
        cam = ctx.camera(image)
        d = cam.ray(px)
        P = np.eye(3) - np.outer(d, d)
        A += P
        b += P @ cam.center
    try:
        X = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        return None, None
    errs = []
    for image, px in sightings:
        cam = ctx.camera(image)
        p = cam.project(X)
        if p is None or cam.depth(X) <= 0:
            return None, None
        errs.append(float(np.linalg.norm(p - np.asarray(px, float))))
    return X, errs


def vet_cluster(ctx, image, cluster, opts):
    """The cluster's usable members when they triangulate cleanly, else ``None``."""
    if cluster["member"]["status"] not in ("reference", "kept"):
        return None
    best: dict[int, dict] = {}
    for m in cluster["members"]:
        if m["image"] < 0 or m["status"] not in ("reference", "kept"):
            continue
        z = m["zncc"] if np.isfinite(m["zncc"]) else 1.0
        if m["image"] not in best or z > best[m["image"]]["_z"]:
            best[m["image"]] = {**m, "_z": z}
    members = list(best.values())
    while len(members) >= 2:
        X, errs = triangulate(ctx, [(m["image"], m["position"]) for m in members])
        if X is not None and max(errs) <= opts["vet_max_reproj_px"]:
            return members
        if X is None or len(members) < 3:
            return None
        worst = int(np.argmax(errs))
        if members[worst]["image"] == image:
            return None
        members.pop(worst)
    return None


def add_vetted(ctx, track, q, image, pixel, opts, info):
    """Carry the pixel through every vetted cluster into the images the track lacks."""
    from sfmtool._sfmtool import bench as B

    near = ctx.clusters_near(image, pixel, opts["vet_radius_px"])
    passed, added = [], []
    seen = {int(o["image"]) for o in track.observations}
    for c in near[: opts["vet_max_clusters"]]:
        members = vet_cluster(ctx, image, c, opts)
        if members is None:
            continue
        passed.append(c["cluster"])
        m = c["member"]
        shape_q = np.asarray(m["shape"], float)
        if abs(np.linalg.det(shape_q)) <= 1e-12:
            continue
        offset = np.asarray(pixel, float) - np.asarray(m["position"], float)
        step = np.linalg.inv(shape_q) @ offset
        for o in members:
            if o["image"] == image or o["image"] in seen:
                continue
            pred = (
                np.asarray(o["position"], float) + np.asarray(o["shape"], float) @ step
            )
            if not in_frame(ctx, o["image"], pred):
                continue
            track, _ = B.add_observation(
                track, int(o["image"]), tuple(map(float, pred)), provenance="sweep"
            )
            seen.add(o["image"])
            added.append(int(o["image"]))
    info["vet"] = {"near": len(near), "passed": passed, "added": added}
    if added:
        track, _ = B.evaluate(track, ctx.edited, ctx.pyramids)
        track, _ = B.apply_thresholds(track)
        track = relax(track, q, opts, info)
        track = anchored_fit(ctx, track, q, pixel, 1)
        track, _ = B.apply_thresholds(track)
    return track


def relax(track, q, opts, info):
    """Turn in the views the thresholds turned out that land where the geometry says."""
    from sfmtool._sfmtool import bench as B

    if not opts["relax"]:
        return track
    in_images = {int(o["image"]) for o in track.observations if o["verdict"] == "in"}
    for i, o in enumerate(track.observations):
        if i == q or o["verdict"] == "in" or int(o["image"]) in in_images:
            continue
        z, s, p = (
            _row(o, k) for k in ("zncc", "seed_shift_px", "projection_offset_px")
        )
        if z is None or s is None or p is None:
            continue
        if (
            z >= opts["relax_min_zncc"]
            and s <= opts["relax_max_shift_px"]
            and p <= opts["relax_max_projection_px"]
        ):
            try:
                track, _ = B.set_verdict(track, i, "in")
            except ValueError:
                continue
            in_images.add(int(o["image"]))
            info.setdefault("relaxed", []).append(int(o["image"]))
    return track


def grow(ctx, track, q, pixel, opts, info, rounds, key):
    from sfmtool._sfmtool import bench as B

    track = relax(track, q, opts, info)
    for _ in range(rounds):
        if track.verdict_counts[0] < 2:
            break
        grown, geo = B.search_geometry(track, q, ctx.edited, ctx.pyramids)
        info.setdefault(key, []).append(geo["added"])
        if not geo["added"]:
            break
        grown, _ = B.evaluate(grown, ctx.edited, ctx.pyramids)
        grown, _ = B.apply_thresholds(grown)
        grown = relax(grown, q, opts, info)
        grown = anchored_fit(ctx, grown, q, pixel, 1)
        track, _ = B.apply_thresholds(grown)
    return track


def _tilted(ctx, track, q, pixel, n):
    from sfmtool._sfmtool import bench as B

    t = shape_track(ctx, track, q, pixel, n)
    t = anchored_fit(ctx, t, q, pixel, 1)
    t, _ = B.apply_thresholds(t)
    return t


def orient(ctx, track, q, image, pixel, opts, ropts, info):
    """The first source of the ``normal`` chain that gives a normal the reading accepts."""
    w, tol = opts["accept_shift_weight"], opts["accept_tolerance"]
    xyz = np.asarray(track.position, float)
    cam = ctx.camera(image)
    to_cam = _unit(cam.center - xyz)
    half = float(np.linalg.norm(track.placement["u_halfvec"]))

    def accept(t, name):
        before, after = reading(track, w), reading(t, w)
        info.setdefault("steps", []).append([name, before, after])
        return after >= before - tol

    for src in opts["normal"].split("+"):
        if src == "nbpos":
            n = renormal.nb_position_normal(ctx, image, pixel, xyz, half, to_cam, ropts)
        elif src == "nb3dpos":
            n = renormal.nb3d_position_normal(ctx, xyz, half, to_cam, ropts)
        elif src in ("photo", "photo_grid", "grid"):
            n = None
            if src != "grid":
                dirs = [
                    _unit(ctx.camera(int(o["image"])).center - xyz)
                    for o in track.observations
                    if o["verdict"] == "in"
                ]
                mean_view = renormal._face(_unit(sum(dirs)), to_cam)
                pn = renormal.photo_normal(ctx, track, mean_view, to_cam, ropts)
                if pn is not None:
                    t = _tilted(ctx, track, q, pixel, pn)
                    if accept(t, "photo"):
                        track, n = t, pn
            if src != "photo" and (
                opts["grid_when"] == "always" or not ctx.points_near(xyz, k=1)
            ):
                ginfo: dict = {}
                gopts = {**cosolve.DEFAULTS, "grid_rounds": opts["grid_rounds"]}
                gn = cosolve.grid_normal(
                    ctx, track, q, pixel, image, gopts, ropts, ginfo
                )
                info["grid"] = ginfo
                if gn is not None:
                    t = _tilted(ctx, track, q, pixel, gn)
                    if accept(t, "grid"):
                        track, n = t, gn
            if n is None:
                continue
            info["source"] = src
            return track
        else:
            raise ValueError(f"unknown normal source {src!r}")
        if n is not None:
            t = _tilted(ctx, track, q, pixel, n)
            # A neighbour's normal is kept whatever the reading says: on a
            # track with few views the right normal can read lower than a wrong
            # one (as renormal's keep="force" found).
            if opts["accept_neighbours"] == "always" or accept(t, src):
                info["source"] = src
                return t
    info["source"] = None
    return track


def rejudge(ctx, track, q, pixel, opts, info):
    """The relaxed views judged again once the normal is settled."""
    from sfmtool._sfmtool import bench as B

    relaxed = set(info.get("relaxed", []))
    if not relaxed:
        return track
    track, _ = B.evaluate(track, ctx.edited, ctx.pyramids)
    dropped = []
    for i, o in enumerate(track.observations):
        if i == q or o["verdict"] != "in" or int(o["image"]) not in relaxed:
            continue
        z = _row(o, "zncc")
        if z is None or z < opts["final_min_zncc"]:
            track, _ = B.set_verdict(track, i, "out")
            dropped.append(int(o["image"]))
    info["rejudged_out"] = dropped
    if dropped and track.verdict_counts[0] >= 2:
        track = anchored_fit(ctx, track, q, pixel, 1)
    return track


def build_track(ctx, image: int, pixel, options: dict | None = None):
    from sfmtool._sfmtool import bench as B

    opts = {**DEFAULTS, **(options or {})}
    ropts = {**renormal.DEFAULTS, "source": "nbpos+nb3dpos+photo"}
    base = renormal.build_track(ctx, image, pixel, {"source": "none"})
    track, q = base.track, base.query_observation
    if track.at_infinity or track.placement is None or q is None:
        return base
    info: dict = {"start_in": int(track.verdict_counts[0])}

    def fallback(reason):
        result = cosolve.build_track(ctx, image, pixel)
        result.diagnostics["staged"] = {**info, "fallback": reason}
        return result

    try:
        if opts["vet_clusters"]:
            track = add_vetted(ctx, track, q, image, pixel, opts, info)
        track = grow(ctx, track, q, pixel, opts, info, opts["grow_rounds"], "grown")
        info["grown_in"] = int(track.verdict_counts[0])
        if track.at_infinity or track.placement is None:
            return fallback("at infinity after growth")
        track = orient(ctx, track, q, image, pixel, opts, ropts, info)
        track = grow(ctx, track, q, pixel, opts, info, opts["regrow_rounds"], "regrown")
        track = rejudge(ctx, track, q, pixel, opts, info)
        track, _ = B.apply_thresholds(track)
        track = clean(ctx, track, q, pixel, ropts, info)
    except ValueError as e:
        return fallback(f"error: {e}")
    if track.at_infinity or not renormal.passes_gates(track, q, pixel, ropts):
        return fallback("gates")
    info["final_in"] = int(track.verdict_counts[0])
    base.diagnostics["staged"] = info
    return TrackAtPixelResult(
        track=track, query_observation=q, diagnostics=base.diagnostics
    )


__all__ = ["build_track", "DEFAULTS", "reading", "triangulate", "vet_cluster"]
