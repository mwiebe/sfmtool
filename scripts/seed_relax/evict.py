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

The rule reads three things the member already holds, and nothing else:

* the observation's **feature radius**, ``0.5 * refine_radius * (|S c0| + |S
  c1|)`` over its stored affine ``S`` -- the reading `fill.source_clusters`
  takes, per OBSERVATION here rather than per cluster, because a footprint is a
  thing in one image: the same cluster is wide on the frame that sees it close
  and narrow on the frame that sees it far;
* its **footprint**, the disk of radius ``2.5 * r / refine_radius`` at its
  pixel position: two and a half of the cluster's own REFINED UNIT SCALES,
  which is the extent the patch cloud states a surfel at, read off the affine
  the matching already refined rather than off a second file;
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

The stage reads the state and the member's own arrays alone; no workspace and
no `.sift` file enter it.  It refuses, and changes nothing, where the state
holds no observations, where the member states no affine shapes, or where the
source file states no refine radius: the footprint is then unstated and a guess
at it would retire evidence on a number nobody measured.
"""

from __future__ import annotations

import numpy as np

from . import fill, pairs, rings, structure
from .fleet_constants import RING_RATIO_P1

#: The half-extent multiplier the footprint is read at, in units of the
#: cluster's own refined unit scale.  Not a tuned number and not a fleet
#: reading: it is the extent the patch cloud states a keypoint's support at, so
#: the disk this rule reads is the disk the surfel occupies.
FOOTPRINT_FACTOR = 2.5

FOOTPRINT_FACTOR_PROVENANCE = {
    "source": "sfmtool_core::patch::PatchExtent::default()",
    "source_file": "crates/sfmtool-core/src/patch/cloud.rs",
    "rule": (
        "`PatchExtent::FeatureSize { factor: 2.5, across: Median }` is the "
        "default sizing policy of the patch cloud: a keypoint carries a patch "
        "of half-extent 2.5 unit scales.  The eviction reads containment at "
        "that same half-extent, stated in the cluster's OWN refined unit "
        "scale `r / refine_radius`, so the footprint it judges is the "
        "footprint the surfel occupies rather than the wider disk the refine "
        "grid measured on, and it is read off the affine the state already "
        "carries rather than off a second file."
    ),
}


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


def footprint(r_px, refine_radius):
    """The observation's drawn footprint in image pixels.

    :data:`FOOTPRINT_FACTOR` of the observation's own refined unit scale, ``r /
    refine_radius``, which is the extent the patch cloud states its support at.
    The refine grid measures on the whole radius ``r``, so a footprint always
    sits inside the disk the pair is read within."""
    return FOOTPRINT_FACTOR * np.asarray(r_px, float) / float(refine_radius)


def covered_by_finer(uv, r_px, refine_radius, slot_i, slot_c, ratio):
    """``(flag, census)``: which observations a finer tracked one covers.

    ``flag`` is one bool per row.  A row is flagged where another row in the
    same image, on another cluster, has its centre inside the row's drawn
    footprint (:func:`footprint`) and a feature radius at least ``ratio`` times
    finer.  The candidates of a footprint are
    :func:`seed_relax.pairs.image_candidates`'s.

    The enumeration is per image and the disk is the coarse row's own
    footprint, so the order rows are visited in cannot change the flag: a row
    is flagged by the existence of a cover, and every cover is measured against
    the state as it stands."""
    uv = np.ascontiguousarray(np.asarray(uv, float))
    r_px = np.asarray(r_px, float)
    slot_i = np.asarray(slot_i, np.int64)
    slot_c = np.asarray(slot_c, np.int64)
    reach = footprint(r_px, refine_radius)
    flag = np.zeros(len(uv), bool)
    n_pairs = 0
    n_octave = 0
    for _img, sel in pairs.image_slices(slot_i):
        if len(sel) < 2:
            continue
        r = r_px[sel]
        c = slot_c[sel]
        rch = reach[sel]
        for big, small, d in pairs.image_candidates(uv[sel, 0], uv[sel, 1], rch):
            inside = (d <= rch[big]) & (r[small] < r[big]) & (c[small] != c[big])
            finer = inside & (r[big] >= ratio * r[small])
            n_pairs += int(inside.sum())
            n_octave += int(finer.sum())
            flag[sel[np.unique(big[finer])]] = True
    census = {
        "n_rows": int(len(uv)),
        "n_pairs_contained": int(n_pairs),
        "n_pairs_finer_band": int(n_octave),
        "n_obs_covered": int(flag.sum()),
    }
    return flag, census


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


def evict_covered(mx, cam, state, refine_radius, floor_px):
    """``(state, census)`` with the covered coarse observations retired.

    Everything the rule reads is the state's own: the rows
    ``structure.state_rows`` names, their stored affines and their cluster ids.
    The state comes back with the covered rows pruned, with the points the
    two-observation rule dropped removed, and adjusted once with the lens
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

    edges = rings.octave_edges(RING_RATIO_P1)
    ratio = octave_ratio(edges)
    r_px = feature_radius(np.asarray(mx.obs_shape, float)[rows], refine_radius)
    flag, pair_census = covered_by_finer(
        np.asarray(mx.obs_uv, float)[rows],
        r_px,
        refine_radius,
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


def evict_stage(mx, cam, state, refine_radius, floor_px, trace=None):
    """The stage as the pipeline runs it: rule, sweep, adjustment."""
    state, census = evict_covered(mx, cam, state, refine_radius, floor_px)
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
