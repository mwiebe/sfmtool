# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Retiring a coarse observation that sits on top of a finer tracked one.

The seed is bootstrapped from a small number of LARGE-radius clusters and every
stage before the fill-in reads only that admission.  Once the fill-in has put
the finer bands back, a coarse feature is no longer the best evidence at the
place it sits: its support spans many image pixels, so it averages over
whatever detail lies under it, and where that detail is at more than one depth
its own triangulated depth is a blend of them.  This stage hands the member
over from the features that carried the bootstrap to the features the fill-in
brought in.

The rule reads three things the member already holds and one the workspace
does:

* the observation's **feature radius**, ``0.5 * refine_radius * (|S c0| + |S
  c1|)`` over its stored affine ``S`` -- the reading `fill.source_clusters`
  takes, per OBSERVATION here rather than per cluster, because a footprint is a
  thing in one image: the same cluster is wide on the frame that sees it close
  and narrow on the frame that sees it far;
* the observation's **keypoint scale** ``sigma``, the column-0 norm of its
  `.sift` affine, which is what the patch cloud sizes a surfel by;
* its **footprint**, the disk of radius ``2.5 sigma`` at its pixel position --
  the support the detector measured, at the extent the patch cloud states it;
* the cluster ids, so a feature never covers itself.

A coarse observation is retired where another observation in the SAME image,
on another cluster the state also holds, has its centre inside that footprint
and a feature radius at least one band finer.  One band is the ratio the band
grid itself is cut on (:func:`seed_relax.rings.octave_edges`), so the rule is a
SCALE rule and never fires on a same-scale neighbour: in a filled-in member
almost every feature has something smaller somewhere inside it, and a rule
without the scale test is a statement about local density rather than about
which evidence is finer.

A point whose surviving observations number fewer than two is dropped, and its
survivors with it, because that is the writer's own rule
(:func:`seed_relax.release.alive_clusters`) and the adjustment cannot constrain
such a point either.  One adjustment with the LENS HELD follows, so the
geometry absorbs the hand-over before the lens is asked.

The stage refuses, and changes nothing, where the workspace's `.sift` files
cannot be read: the footprint is then unknown and a guess at it would retire
evidence on a number nobody measured.
"""

from __future__ import annotations

from pathlib import Path, PurePosixPath

import numpy as np

from . import fill, rings, structure
from .fleet_constants import RING_RATIO_P1

#: The half-extent multiplier the footprint is read at.  Not a tuned number and
#: not a fleet reading: it is the extent the patch cloud states a keypoint's
#: support at, so the disk this rule reads is the disk the surfel occupies.
FOOTPRINT_FACTOR = 2.5

FOOTPRINT_FACTOR_PROVENANCE = {
    "source": "sfmtool_core::patch::PatchExtent::default()",
    "source_file": "crates/sfmtool-core/src/patch/cloud.rs",
    "rule": (
        "`PatchExtent::FeatureSize { factor: 2.5, across: Median }` is the "
        "default sizing policy of the patch cloud: a keypoint of scale sigma "
        "carries a patch of half-extent 2.5 sigma, sigma being the column-0 "
        "norm of the `.sift` affine shape.  The eviction reads containment at "
        "that same half-extent, so the footprint it judges is the footprint "
        "the surfel occupies rather than the wider disk the refine grid "
        "measured on."
    ),
}

#: How many candidate (large, small) row pairs one pass of the containment
#: expansion holds at a time.  A memory bound on a batched computation, not a
#: threshold: the pairs are the same pairs and the flag is the same flag at any
#: value of it.
_PAIR_CHUNK = 1 << 22


def evict_on(env):
    """Whether the eviction stage runs (``SFMTOOL_RELAX_EVICT``)."""
    return (env.get("SFMTOOL_RELAX_EVICT", "1") or "1").strip() != "0"


def octave_ratio(edges):
    """The radius ratio one band spans, read off the band grid itself.

    ``edges`` is :func:`seed_relax.rings.octave_edges`'s output, whose top edge
    is infinite; the first two finite edges state the ratio every band below
    them is cut on, so "one band finer" costs no constant of its own."""
    finite = [float(e) for e in edges if np.isfinite(e) and e > 0.0]
    if len(finite) < 2:
        raise ValueError(f"band grid states no ratio: {edges!r}")
    return finite[0] / finite[1]


def feature_radius(obs_shape, refine_radius):
    """The observation's feature radius in image pixels.

    `fill.source_clusters`'s reading, per observation: the refine radius times
    the mean of the stored affine's two column norms."""
    sh = np.asarray(obs_shape, float)
    return (
        0.5
        * float(refine_radius)
        * (np.linalg.norm(sh[:, :, 0], axis=1) + np.linalg.norm(sh[:, :, 1], axis=1))
    )


def sift_path(workspace, feature_prefix_dir, name):
    """The `.sift` file beside an image, in the workspace's own layout.

    ``{workspace}/{image parent}/{feature_prefix_dir}/{image name}.sift``, the
    path :meth:`SfmrReconstruction::sift_path_for_image` states, so a rig
    capture whose frames sit under one directory per camera resolves each
    camera's own feature directory."""
    rel = PurePosixPath(str(name).replace("\\", "/"))
    return Path(workspace) / rel.parent / str(feature_prefix_dir) / f"{rel.name}.sift"


def keypoint_scales(names, obs_i, obs_f, workspace, feature_prefix_dir):
    """``sigma`` per observation: the column-0 norm of its `.sift` affine.

    The same reading the patch cloud takes.  Each image's file is read once,
    capped at the highest feature index the rows name.  Returns ``None`` where
    the workspace does not state a feature directory; raises where a file the
    rows name cannot be read."""
    if not workspace or not feature_prefix_dir:
        return None
    from sfmtool._sfmtool.io import read_sift_partial

    obs_i = np.asarray(obs_i, np.int64)
    obs_f = np.asarray(obs_f, np.int64)
    sig = np.full(len(obs_i), np.nan)
    for img in np.unique(obs_i):
        sel = obs_i == img
        need = int(obs_f[sel].max()) + 1
        path = sift_path(workspace, feature_prefix_dir, names[int(img)])
        data = read_sift_partial(str(path), need)
        aff = np.asarray(data["affine_shapes"], float)
        sig[sel] = np.hypot(aff[obs_f[sel], 0, 0], aff[obs_f[sel], 1, 0])
    return sig


def covered_by_finer(uv, r_px, sigma, slot_i, slot_c, ratio):
    """``(flag, census)``: which observations a finer tracked one covers.

    ``flag`` is one bool per row.  A row is flagged where another row in the
    same image, on another cluster, has its centre inside the row's drawn
    footprint (``2.5 sigma``, and inside the feature radius, which is the
    support the pair is read within) and a feature radius at least ``ratio``
    times finer.

    The enumeration is per image and the disk is the smaller of the two radii,
    so the pair set is the conjunction of both containments and the order rows
    are visited in cannot change the flag: a row is flagged by the existence of
    a cover, and every cover is measured against the state as it stands."""
    uv = np.ascontiguousarray(np.asarray(uv, float))
    r_px = np.asarray(r_px, float)
    sigma = np.asarray(sigma, float)
    slot_i = np.asarray(slot_i, np.int64)
    slot_c = np.asarray(slot_c, np.int64)
    reach = np.minimum(r_px, FOOTPRINT_FACTOR * sigma)
    flag = np.zeros(len(uv), bool)
    n_pairs = 0
    n_octave = 0
    for img in np.unique(slot_i):
        sel = np.nonzero(slot_i == img)[0]
        if len(sel) < 2:
            continue
        x = np.ascontiguousarray(uv[sel, 0])
        y = np.ascontiguousarray(uv[sel, 1])
        r = r_px[sel]
        c = slot_c[sel]
        rch = reach[sel]
        # The rows of this image ordered by x, so the candidates of a disk are
        # a contiguous run: a disk of radius `rch` cannot reach past a column
        # further than `rch` away.
        xo = np.argsort(x, kind="stable")
        xs = x[xo]
        lo = np.searchsorted(xs, x - rch, side="left")
        hi = np.searchsorted(xs, x + rch, side="right")
        cnt = (hi - lo).astype(np.int64)
        cnt[~np.isfinite(rch)] = 0
        for a, b in _chunks(cnt):
            got = _pairs(x, y, r, c, rch, xo, lo, cnt, a, b, ratio)
            n_pairs += int(got[0])
            n_octave += int(got[1])
            flag[sel[got[2]]] = True
    census = {
        "n_rows": int(len(uv)),
        "n_pairs_contained": int(n_pairs),
        "n_pairs_finer_band": int(n_octave),
        "n_obs_covered": int(flag.sum()),
    }
    return flag, census


def _chunks(cnt):
    """``(start, stop)`` row ranges whose candidate pairs fit one pass."""
    ends = np.cumsum(cnt)
    a = 0
    while a < len(cnt):
        base = ends[a - 1] if a else 0
        b = int(np.searchsorted(ends, base + _PAIR_CHUNK, side="right"))
        b = max(b, a + 1)
        yield a, min(b, len(cnt))
        a = b


def _pairs(x, y, r, c, rch, xo, lo, cnt, a, b, ratio):
    """One chunk of the expansion: ``(n_pairs, n_finer, flagged rows)``."""
    take = cnt[a:b]
    total = int(take.sum())
    if not total:
        return 0, 0, np.zeros(0, np.int64)
    big = np.repeat(np.arange(a, b, dtype=np.int64), take)
    offs = np.arange(total, dtype=np.int64) - np.repeat(np.cumsum(take) - take, take)
    small = xo[np.repeat(lo[a:b].astype(np.int64), take) + offs]
    dx = x[small] - x[big]
    dy = y[small] - y[big]
    d = np.sqrt(dx * dx + dy * dy)
    inside = (d <= rch[big]) & (r[small] < r[big]) & (c[small] != c[big])
    finer = inside & (r[big] >= ratio * r[small])
    return int(inside.sum()), int(finer.sum()), np.unique(big[finer])


def two_observation_sweep(slot_c, keep_obs, n_points):
    """``(keep_obs, keep_point)`` with the writer's own two-observation rule.

    A point whose surviving observations number fewer than two is dropped and
    its survivors go with it, so nothing downstream reads a row whose point the
    writer would not keep."""
    slot_c = np.asarray(slot_c, np.int64)
    keep_obs = np.asarray(keep_obs, bool)
    keep_pt = np.bincount(slot_c[keep_obs], minlength=int(n_points)) >= 2
    return keep_obs & keep_pt[slot_c], keep_pt


def band_census(ring, keep_obs, slot_c, keep_pt):
    """Per band: how many observations and points went in and came out."""
    ring = np.asarray(ring, np.int64)
    if not len(ring):
        return []
    # A cluster's band is its WIDEST observation's, which is the reading the
    # fill-in banded it on.
    cl_ring = np.full(len(keep_pt), np.iinfo(np.int64).max, np.int64)
    np.minimum.at(cl_ring, np.asarray(slot_c, np.int64), ring)
    out = []
    for r in range(int(ring.max()) + 1):
        rows = ring == r
        pts = cl_ring == r
        if not rows.any() and not pts.any():
            continue
        out.append(
            {
                "ring": int(r),
                "obs": int(rows.sum()),
                "obs_kept": int((rows & keep_obs).sum()),
                "points": int(pts.sum()),
                "points_kept": int((pts & keep_pt).sum()),
            }
        )
    return out


def evict_covered(mx, cam, state, sigma, refine_radius, floor_px):
    """``(state, census)`` with the covered coarse observations retired.

    ``sigma`` is one keypoint scale per row of ``structure.state_rows``, in
    that order.  The state comes back with those rows pruned, with the points
    the two-observation rule dropped removed, and adjusted once with the lens
    held."""
    rows, _slot_i, slot_c = structure.state_rows(mx, state)
    n_pts = len(np.asarray(state["clusters"]))
    census = {"footprint_factor": FOOTPRINT_FACTOR, "n_rows": int(len(rows))}
    if not len(rows):
        census["refused"] = "the state holds no observations"
        return state, census
    if mx.obs_shape is None:
        census["refused"] = "the member states no affine shapes"
        return state, census
    if not refine_radius:
        census["refused"] = "the source file states no refine radius"
        return state, census
    sigma = np.asarray(sigma, float)
    if len(sigma) != len(rows) or not np.isfinite(sigma).all():
        census["refused"] = "a keypoint scale could not be read"
        return state, census

    edges = rings.octave_edges(RING_RATIO_P1)
    ratio = octave_ratio(edges)
    r_px = feature_radius(np.asarray(mx.obs_shape, float)[rows], refine_radius)
    flag, pair_census = covered_by_finer(
        np.asarray(mx.obs_uv, float)[rows],
        r_px,
        sigma,
        np.asarray(mx.obs_i, np.int64)[rows],
        np.asarray(mx.obs_c, np.int64)[rows],
        ratio,
    )
    census.update(pair_census)
    census["band_ratio"] = float(ratio)
    keep_obs, keep_pt = two_observation_sweep(slot_c, ~flag, n_pts)

    n_obs_pt = np.bincount(slot_c, minlength=n_pts)
    n_ev_pt = np.bincount(slot_c[flag], minlength=n_pts)
    dropped = ~keep_pt
    census.update(
        {
            "n_obs_evicted": int((~keep_obs).sum()),
            "n_obs_kept": int(keep_obs.sum()),
            "n_points": int(n_pts),
            "n_points_kept": int(keep_pt.sum()),
            "n_points_dropped": int(dropped.sum()),
            "n_points_dropped_all_covered": int(
                (dropped & (n_ev_pt >= n_obs_pt)).sum()
            ),
            "n_points_dropped_by_two_obs": int((dropped & (n_ev_pt < n_obs_pt)).sum()),
        }
    )
    if floor_px:
        census["bands"] = band_census(
            rings.assign_rings(r_px, float(floor_px), edges),
            keep_obs,
            slot_c,
            keep_pt,
        )
    if keep_obs.all():
        census["adjusted"] = False
        return state, census

    state = structure.with_pruned(state, rows[~keep_obs])
    if not keep_pt.all():
        state = dict(state)
        state["clusters"] = np.asarray(state["clusters"], np.int64)[keep_pt]
        state["points"] = np.asarray(state["points"], float)[keep_pt]
        state["at_inf"] = np.asarray(state["at_inf"], bool)[keep_pt]
    state, barec = fill.adjust_held(mx, cam, state)
    census["adjusted"] = True
    census["ba_n_obs"] = barec["n_obs"]
    census["ba_reproj_med_px"] = barec["reproj_med_px"]
    census["ba_reproj_p90_px"] = barec["reproj_p90_px"]
    return state, census


def source_scales(mx, state, source, workspace):
    """``(sigma, refused)`` for the state's own rows, off the workspace.

    The feature directory is the `.matches` file's own record of it, so the
    scales come from the files the capture was matched on rather than from a
    directory guessed at."""
    if not workspace:
        return None, "no workspace given"
    meta = getattr(source, "metadata", None) or {}
    prefix = ((meta.get("workspace") or {}).get("contents") or {}).get(
        "feature_prefix_dir"
    )
    if not prefix:
        return None, "the source file names no feature directory"
    rows, _slot_i, _slot_c = structure.state_rows(mx, state)
    try:
        sig = keypoint_scales(
            list(mx.names),
            np.asarray(mx.obs_i, np.int64)[rows],
            np.asarray(mx.obs_f, np.int64)[rows],
            workspace,
            prefix,
        )
    except Exception as exc:  # noqa: BLE001 -- an unreadable file is a refusal
        return None, f"{type(exc).__name__}: {exc}"
    return sig, None


def evict_stage(mx, cam, state, source, workspace, refine_radius, floor_px, trace=None):
    """The stage as the pipeline runs it: scales, rule, sweep, adjustment."""
    sig, why = source_scales(mx, state, source, workspace)
    if sig is None:
        why = why or "the keypoint scales could not be read"
        census = {"footprint_factor": FOOTPRINT_FACTOR, "refused": why}
        if trace is not None:
            trace(f"    evict: refused ({why})")
        return state, census
    state, census = evict_covered(mx, cam, state, sig, refine_radius, floor_px)
    if trace is not None:
        if census.get("refused"):
            trace(f"    evict: refused ({census['refused']})")
        else:
            trace(
                f"    evict: {census['n_obs_evicted']} observations retired of "
                f"{census['n_rows']}, {census['n_points_kept']} points kept of "
                f"{census['n_points']}, reproj med "
                f"{census.get('ba_reproj_med_px')}"
            )
    return state, census
