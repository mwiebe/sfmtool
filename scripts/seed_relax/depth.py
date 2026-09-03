# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""A depth is stated only where one was measured.

The re-estimation before this stage asks whether a point's rays cross at all:
a point whose widest ray pair does not clear the member's angular floor has no
measured depth and becomes a bearing.  That floor is a statement about
DISTINGUISHABILITY, that two rays are further apart than the pixel bar, and
not about how well the crossing is placed.  A pair barely past the floor crosses at
a depth the same pixel bar moves by most of the distance to the point, so the
member ships a position it cannot tell from any other position along the ray.

This stage reads that error and demotes the points whose depth it swallows.
For every finite point at the final geometry:

* ``theta_max`` is the widest angle between any pair of its world rays, the
  same statistic the floor rule reads (the kernel reads it per point but does
  not state it, so it is read again here from the same rays at the same
  poses and the same lens);
* ``eps`` is the point's own median reprojection residual in pixels, floored at
  the member's median over its finite points, so a lucky two-observation track
  is not trusted below what the member as a whole achieves;
* ``u = eps / (f_eq * sin theta_max)`` is the RELATIVE depth uncertainty, the
  pixel error carried through the equivalent focal and divided by the parallax
  that resolves it;
* ``du = u * |Z|`` is that error at the point's own depth ``Z``, its median
  depth over the cameras that see it.

The bound is the member's own local support radius: ``r12``, the median
distance to the :data:`K_SUPPORT`-th nearest neighbour among the CONFIDENT HALF
of the finite points, those whose ``u`` is at or below its median.  A point
whose depth error is larger than the radius of the neighbourhood its own best
points define cannot be placed inside that neighbourhood: whatever surface it
samples, it is not a sample of that surface at a known position.

A point over the bound becomes a BEARING: ``w = 0``, direction the normalised
mean of its world rays.  Nothing is deleted: every point and every
observation survives the stage, and only the depth claim is withdrawn, so a
demoted point still carries its direction, still constrains rotation, and is
still there to be re-read by anything downstream that measures it again.

Everything the rule reads is the member's own: its residuals, its rays, its
equivalent focal and the spacing of its own confident points.  No reference, no
capture-level constant and no absolute length enters it.
"""

from __future__ import annotations

import numpy as np

from . import structure

#: The neighbourhood size the support radius is read at.  Not a bound of this
#: rule's own: it is the neighbour count the patch cloud and the ARS normal
#: estimator already read a point's local surface in, so ``r12`` is the radius
#: of the patch of surface a point belongs to as the rest of the toolkit
#: already measures it.
K_SUPPORT = 12

K_SUPPORT_PROVENANCE = {
    "source": "sfmtool.reconstruction.SfmrReconstruction.to_embedded_patches",
    "source_arg": "k_neighbors=12",
    "also": "sfmtool_core::patch::normals (ARS adjacency neighbourhood)",
    "rule": (
        "Twelve is the neighbourhood the embed and the ARS normal estimator "
        "read a point's local surface in.  The support radius this stage "
        "compares a depth error against is the radius of that same "
        "neighbourhood, measured on the points whose depth is best determined, "
        "so the length costs no constant of its own and means the same thing "
        "here as it does wherever else a neighbourhood is read."
    ),
}

#: The angle bins the census reports the demotion over, in degrees.
ANGLE_EDGES = (0.0, 1.0, 2.0, 4.0, 8.0, np.inf)
ANGLE_NAMES = ("under 1 deg", "1 to 2 deg", "2 to 4 deg", "4 to 8 deg", "over 8 deg")

#: How many ray pairs one pass of the widest-angle scan holds at a time.  A
#: memory bound and not a threshold: the angles are the same angles at any
#: value of it.
PAIR_CHUNK = 1 << 22


# ------------------------------------------------------------------ per point


def groups(slot_c, n_points):
    """``(order, counts, starts)``: the rows of each point, contiguous.

    ``order`` sorts the rows by the point they name, stably, so a point's rows
    are ``order[starts[c] : starts[c] + counts[c]]`` and their order inside the
    block is the order the state holds them in."""
    slot_c = np.asarray(slot_c, np.int64)
    order = np.argsort(slot_c, kind="stable")
    counts = np.bincount(slot_c, minlength=int(n_points)).astype(np.int64)
    starts = np.concatenate([[0], np.cumsum(counts)])[:-1].astype(np.int64)
    return order, counts, starts


def _blocks(points, count, chunk=PAIR_CHUNK):
    """Chunk a point list so one pass holds at most ``chunk`` ray pairs."""
    step = max(1, int(chunk // (int(count) * int(count) + 1)))
    for lo in range(0, len(points), step):
        yield points[lo : lo + step]


def widest_ray_angles(rays, order, counts, starts):
    """The widest angle between any pair of a point's rays, in radians.

    Zero where a point has fewer than two observations, which is the angle a
    single ray subtends with itself.  ``rays`` are the unit world rays of the
    state's rows, in row order."""
    rays = np.asarray(rays, float)
    out = np.zeros(len(counts))
    for c in np.unique(counts):
        if c < 2:
            continue
        pts = np.nonzero(counts == c)[0]
        for blk in _blocks(pts, int(c)):
            idx = starts[blk][:, None] + np.arange(int(c))[None, :]
            r = rays[order[idx]]
            cos = np.einsum("mik,mjk->mij", r, r)
            out[blk] = np.arccos(np.clip(cos.min(axis=(1, 2)), -1.0, 1.0))
    return out


def point_medians(values, order, counts, starts):
    """The median of a per-observation quantity over each point's own rows.

    ``NaN`` where a point has no observation at all."""
    v = np.asarray(values, float)[order]
    out = np.full(len(counts), np.nan)
    for c in np.unique(counts):
        if c < 1:
            continue
        pts = np.nonzero(counts == c)[0]
        for blk in _blocks(pts, int(c)):
            idx = starts[blk][:, None] + np.arange(int(c))[None, :]
            out[blk] = np.median(v[idx], axis=1)
    return out


def ray_depths(points, centres, rays, slot_c, slot_i):
    """Each observation's depth of its own point along its own ray."""
    p = np.asarray(points, float)[np.asarray(slot_c, np.int64)]
    c = np.asarray(centres, float)[np.asarray(slot_i, np.int64)]
    return np.einsum("nj,nj->n", p - c, np.asarray(rays, float))


def mean_bearings(rays, order, counts, starts, which):
    """The unit mean of each named point's own world rays.

    One row per point of ``which``, in that array's order; a point whose rays
    cancel keeps whatever direction the mean states and is normalised by the
    same rule as the rest, which the caller's own guard on the norm covers."""
    rays = np.asarray(rays, float)
    which = np.asarray(which, np.int64)
    out = np.zeros((len(which), 3))
    if not len(which):
        return out
    slot = np.full(len(counts), -1, np.int64)
    slot[which] = np.arange(len(which), dtype=np.int64)
    for c in np.unique(counts[which]):
        if c < 1:
            continue
        pts = which[counts[which] == c]
        for blk in _blocks(pts, int(c)):
            idx = starts[blk][:, None] + np.arange(int(c))[None, :]
            d = rays[order[idx]].mean(axis=1)
            d = d / np.maximum(np.linalg.norm(d, axis=1, keepdims=True), 1e-300)
            out[slot[blk]] = d
    return out


# ------------------------------------------------------------------ the bound


def support_radius(positions, k=K_SUPPORT):
    """The median distance to the ``k``-th nearest neighbour of a point set.

    ``None`` where the set holds fewer than ``k + 1`` points: the neighbourhood
    is then not stated, and a radius read off a smaller one would be a reading
    of how few points there are rather than of how far apart they sit."""
    from sfmtool._sfmtool.spatial import KdTree3d

    p = np.ascontiguousarray(np.asarray(positions, np.float64))
    if len(p) < int(k) + 1:
        return None
    idx = np.asarray(KdTree3d(p).nearest_k(p, int(k) + 1), np.int64)
    far = idx[:, -1]
    d = np.linalg.norm(p[far] - p, axis=1)
    return float(np.median(d))


def uncertainties(theta, eps, f_eq):
    """``u = eps / (f_eq * sin theta_max)``, the relative depth uncertainty.

    ``NaN`` where the angle states no parallax, which is a point the floor rule
    has already made a bearing."""
    s = np.sin(np.clip(np.asarray(theta, float), 0.0, 0.5 * np.pi))
    s = np.where(s > 0.0, s, np.nan)
    return np.asarray(eps, float) / (float(f_eq) * s)


def angle_bins(theta, finite, demoted):
    """Per angle bin: how many finite points went in and how many were demoted."""
    deg = np.degrees(np.asarray(theta, float))
    rows = []
    for name, lo, hi in zip(ANGLE_NAMES, ANGLE_EDGES[:-1], ANGLE_EDGES[1:]):
        b = finite & (deg >= lo) & (deg < hi)
        n = int(b.sum())
        d = int((b & demoted).sum())
        rows.append(
            {
                "bin": name,
                "n": n,
                "demoted": d,
                "frac": (d / n) if n else None,
            }
        )
    return rows


# ------------------------------------------------------------------ the stage


def demote_uncertain(mx, cam, state, f_eq, trace=None):
    """``(state, census)`` with the points whose depth error swallows them.

    The state comes back with those points at ``w = 0``, pointing along the
    unit mean of their own world rays.  Every point and every observation the
    state held is still in it; only the depth claim is withdrawn.

    The stage REFUSES, and changes nothing, where the state holds no
    observations, where it states no finite point, where the member states no
    equivalent focal, or where the confident half is smaller than the
    neighbourhood the support radius is read in: the bound is then unstated,
    and a guess at it would withdraw a measurement on a number nobody made."""
    census = {"k_support": int(K_SUPPORT)}
    rows, slot_i, slot_c = structure.state_rows(mx, state)
    at_inf = np.asarray(state["at_inf"], bool)
    points = np.asarray(state["points"], float)
    n_pts = int(len(at_inf))
    census.update(
        {
            "n_rows": int(len(rows)),
            "n_points": n_pts,
            "n_finite_before": int((~at_inf).sum()),
            "n_bearings_before": int(at_inf.sum()),
            "n_demoted": 0,
            "n_finite_after": int((~at_inf).sum()),
            "n_bearings_after": int(at_inf.sum()),
        }
    )
    if not len(rows):
        census["refused"] = "the state holds no observations"
        return state, census
    if not (~at_inf).any():
        census["refused"] = "the state states no finite point"
        return state, census
    if not (np.isfinite(f_eq) and float(f_eq) > 0.0):
        census["refused"] = "the member states no equivalent focal"
        return state, census

    rot, cen = structure.centres_of(state)
    uv = np.asarray(mx.obs_uv, float)[rows]
    rays = structure.world_rays(cam, rot, uv, slot_i)
    order, counts, starts = groups(slot_c, n_pts)
    theta = widest_ray_angles(rays, order, counts, starts)

    resid = structure.reprojection(cam, rot, cen, points, at_inf, uv, slot_i, slot_c)
    fin_row = ~at_inf[slot_c]
    finite_res = resid[fin_row & np.isfinite(resid)]
    if not len(finite_res):
        census["refused"] = "the state states no finite residual"
        return state, census
    eps_member = float(np.median(finite_res))
    eps = np.fmax(point_medians(resid, order, counts, starts), eps_member)

    depth = point_medians(
        ray_depths(points, cen, rays, slot_c, slot_i), order, counts, starts
    )
    u = uncertainties(theta, eps, f_eq)
    finite = (~at_inf) & np.isfinite(u) & np.isfinite(depth)
    census.update(
        {
            "f_eq_px": float(f_eq),
            "eps_member_px": eps_member,
            "n_finite_read": int(finite.sum()),
            "theta_deg_p10": _pct(np.degrees(theta[finite]), 10),
            "theta_deg_p50": _pct(np.degrees(theta[finite]), 50),
            "theta_deg_p90": _pct(np.degrees(theta[finite]), 90),
            "u_p10": _pct(u[finite], 10),
            "u_p50": _pct(u[finite], 50),
            "u_p90": _pct(u[finite], 90),
        }
    )
    if not finite.any():
        census["refused"] = "no finite point states a parallax"
        return state, census

    confident = finite & (u <= float(np.median(u[finite])))
    census["n_confident"] = int(confident.sum())
    r12 = support_radius(points[confident])
    if r12 is None:
        census["refused"] = (
            f"the confident half holds fewer than {K_SUPPORT + 1} points"
        )
        return state, census
    z_conf = float(np.median(np.abs(depth[confident])))
    census.update(
        {
            "support_r12": float(r12),
            "depth_med_confident": z_conf,
            "support_scalar": (float(r12) / z_conf) if z_conf > 0.0 else None,
        }
    )

    demoted = finite & (u * np.abs(depth) > float(r12))
    census["bins"] = angle_bins(theta, finite, demoted)
    census["n_demoted"] = int(demoted.sum())
    census["n_finite_after"] = int((~at_inf).sum() - demoted.sum())
    census["n_bearings_after"] = int(at_inf.sum() + demoted.sum())
    if trace is not None:
        trace(
            f"    depth: {census['n_demoted']} of {census['n_finite_before']} "
            f"finite points demoted past r12 {r12:.6g} "
            f"({census['support_scalar']:.4g} of the confident median depth), "
            f"{census['n_finite_after']} finite left"
        )
    if not demoted.any():
        return state, census

    which = np.nonzero(demoted)[0]
    out = dict(state)
    pts = points.copy()
    pts[which] = mean_bearings(rays, order, counts, starts, which)
    inf = at_inf.copy()
    inf[which] = True
    out["points"] = pts
    out["at_inf"] = inf
    return out, census


def depth_stage(mx, cam, state, f_eq, trace=None):
    """The stage as the pipeline runs it: read the error, withdraw the depth."""
    state, census = demote_uncertain(mx, cam, state, f_eq, trace)
    if trace is not None and census.get("refused"):
        trace(f"    depth: refused ({census['refused']})")
    return state, census


def _pct(v, q):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(np.percentile(v, q)) if len(v) else None
