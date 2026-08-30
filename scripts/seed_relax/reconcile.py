# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Several tracked points resting on one measurement, reconciled into one.

SIFT emits one stored feature per dominant orientation at a detection, so one
corner at one scale enters the feature file two, three or four times with the
same pixel and the same affine scale.  The matching gives every descriptor to
at most one cluster, so those rows are matched into DIFFERENT clusters, and the
relaxation then tracks two or three points that rest on the same physical
measurement.  Nothing downstream can tell them apart: each looks like an
independent track, so a rule that thins observations sees several fragile
points where the evidence states one.

This stage reads the state's own rows, groups the ones that measure one place
in one image, and asks the arithmetic what to do with each group:

* a KNOT is a connected component over "these two points share a measurement";
* its UNION TRACK is one ray per distinct detection over every row of every
  member, solved at the member's own poses and lens;
* the knot is merged into one point, thinned by the one member the resolve
  votes against, or left exactly as it stood.

Every bar is the member's own.  The estimation floors at the member's angular
consensus tolerance, the residual test is PAIRED against the same rows' fit
under the separate points, and the position test is that same tolerance carried
to depth through each member's own widest ray pair.  Nothing absolute is
chosen.

The stage sits between the fill-in and the hand-over, so both the hand-over's
two-observation sweep and the re-estimation's cheirality-minority prune read
merged tracks: a track that carries fifteen detections is not thinned into
three unchecked pairs, and a minority vote is taken over more evidence.
"""

from __future__ import annotations

import math

import numpy as np

from . import evict, pairs, structure

#: Radius agreement two detections need before they are read as one place: the
#: difference over the larger of the two.  A tenth is the survey's own bar,
#: carried here because a wider one starts joining a feature to the coarser one
#: it sits inside, which is the hand-over's relation and not this one.
NEAR_RADIUS_FRAC = 0.10

#: The positional tolerance, as a fraction of the SMALLER row's refined unit
#: scale ``r / refine_radius``.  A fixed pixel bar is not scale aware: the same
#: bar is a fifth of a coarse feature's unit scale and several times a fine
#: one's, so it joins fine detections a whole feature apart while refusing
#: coarse ones that are the same corner.  Stated in unit scales the relation
#: reads the same in every band.
NEAR_UNIT_FRAC = 0.60

NEAR_UNIT_FRAC_PROVENANCE = {
    "fleet": "evo-survey-20260823 / relax-20260827 shared-observation survey",
    "source_table": "shared/prod/tolerance_sweep.json",
    "members": ("20250906_211742965 h14", "KerryPark360 h12"),
    "rule": (
        "The fraction whose pair set best reproduces the survey's fixed 1 px "
        "set on the two deep members, by Jaccard index over unordered row "
        "pairs under the same radius agreement, swept from 0.05 to 3.0 in "
        "steps of 0.05.  Both members peak at 0.60 independently (h14 0.947 "
        "over 17811 reference pairs, Kerry 0.964 over 3261), pooled 0.950; "
        "the optimum is flat, holding above 0.94 pooled from 0.50 to 0.90 and "
        "falling to 0.88 at 0.30 and to 0.63 at 3.0."
    ),
}

#: `estimate_points` verdicts that state a depth.  A rescued cluster is finite;
#: the kernel keeps the two buckets apart so its verdicts still partition.
FINITE_VERDICTS = (0, 6)
#: The verdict that states a direction and no depth.
THIN_VERDICT = 2


def reconcile_on(env):
    """Whether the reconciliation stage runs (``SFMTOOL_RELAX_RECONCILE``)."""
    return (env.get("SFMTOOL_RELAX_RECONCILE", "1") or "1").strip() != "0"


def _ev():
    """The evaluation battery's member class, imported on use."""
    import seed_candidate_eval as EV

    return EV


# ----------------------------------------------------------------- the relation


def detection_groups(slot_i, uv):
    """One group id per ``(image, exact stored pixel)``: detection identity.

    The pixels are compared by their stored bits, so two rows are one detection
    exactly when the feature file put them at the same place.  That is the
    degenerate case of the near relation below, and it needs no tolerance."""
    b = np.ascontiguousarray(np.asarray(uv, np.float64)).view(np.int64).reshape(-1, 2)
    keys = np.ascontiguousarray(
        np.stack([np.asarray(slot_i, np.int64), b[:, 0], b[:, 1]], axis=1)
    )
    view = keys.view(np.dtype((np.void, keys.dtype.itemsize * keys.shape[1]))).ravel()
    _u, gid = np.unique(view, return_inverse=True)
    return np.asarray(gid, np.int64).ravel()


def near_rows(
    slot_i,
    slot_c,
    uv,
    r_px,
    refine_radius,
    frac=NEAR_UNIT_FRAC,
    r_frac=NEAR_RADIUS_FRAC,
):
    """``(a, b)`` row pairs that measure one place in one image.

    Two rows of DIFFERENT clusters in one image are one place when their
    centres are within ``frac`` of the smaller row's refined unit scale and
    their feature radii agree to ``r_frac`` of the larger.  Rows at the exact
    same pixel are not returned: :func:`detection_groups` already states them,
    and it does so without a tolerance.

    Each unordered pair once, in row order."""
    uv = np.asarray(uv, float)
    slot_c = np.asarray(slot_c, np.int64)
    r_px = np.asarray(r_px, float)
    unit = r_px / float(refine_radius)
    A, B = [], []
    for _img, sel in pairs.image_slices(slot_i):
        if len(sel) < 2:
            continue
        reach = float(frac) * unit[sel]
        c = slot_c[sel]
        r = r_px[sel]
        for big, small, d in pairs.image_candidates(uv[sel, 0], uv[sel, 1], reach):
            keep = (
                (d <= reach[big])
                & (d <= reach[small])
                & (d > 0.0)
                & (c[big] != c[small])
            )
            if not keep.any():
                continue
            ra, rb = r[big[keep]], r[small[keep]]
            agree = np.abs(ra - rb) <= float(r_frac) * np.maximum(ra, rb)
            if not agree.any():
                continue
            i = sel[big[keep][agree]]
            j = sel[small[keep][agree]]
            A.append(np.minimum(i, j))
            B.append(np.maximum(i, j))
    if not A:
        z = np.zeros(0, np.int64)
        return z, z
    a = np.concatenate(A).astype(np.int64)
    b = np.concatenate(B).astype(np.int64)
    key = np.unique(a * (1 << 32) + b)
    return (key >> 32).astype(np.int64), (key & ((1 << 32) - 1)).astype(np.int64)


def _settled(parent):
    """Union-find roots, squashed by pointer doubling to a fixed point."""
    while True:
        nxt = parent[parent]
        if np.array_equal(nxt, parent):
            return parent
        parent = nxt


def knot_components(slot_c, gid, pair_a, pair_b, n_points):
    """``(knot per point, knot count, shared row flag)`` over the relation.

    A point is in a knot when one of its rows shares a detection with a row of
    another point, or is a near pair's partner.  Every union carries the lower
    point id upward, and knots are then numbered in point order, so the
    numbering is a function of the state and not of the order the pairs arrived
    in.  ``-1`` marks a point outside every knot."""
    slot_c = np.asarray(slot_c, np.int64)
    gid = np.asarray(gid, np.int64)
    parent = np.arange(int(n_points), dtype=np.int64)
    shared = np.zeros(len(slot_c), bool)

    edges = []
    order = np.argsort(gid, kind="stable")
    g = gid[order]
    cuts = np.flatnonzero(np.diff(g)) + 1
    for lo, hi in zip(np.concatenate(([0], cuts)), np.concatenate((cuts, [len(g)]))):
        if hi - lo < 2:
            continue
        rr = order[lo:hi]
        pts = np.unique(slot_c[rr])
        if len(pts) < 2:
            continue
        shared[rr] = True
        edges.append(np.stack([np.full(len(pts) - 1, pts[0]), pts[1:]], axis=1))
    if len(pair_a):
        shared[pair_a] = True
        shared[pair_b] = True
        edges.append(np.stack([slot_c[pair_a], slot_c[pair_b]], axis=1))

    joins = np.concatenate(edges) if edges else np.zeros((0, 2), np.int64)
    for a, b in joins:
        ra, rb = int(a), int(b)
        while parent[ra] != ra:
            ra = int(parent[ra])
        while parent[rb] != rb:
            rb = int(parent[rb])
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    parent = _settled(parent)

    involved = np.zeros(int(n_points), bool)
    involved[slot_c[shared]] = True
    knot = np.full(int(n_points), -1, np.int64)
    seen = {}
    for p in np.nonzero(involved)[0]:
        r = int(parent[p])
        if r not in seen:
            seen[r] = len(seen)
        knot[p] = seen[r]
    return knot, len(seen), shared


# -------------------------------------------------------------- the union solve


def representatives(gid, slot_c, n_obs_pt, rows_sel=None):
    """``(group, row)``: the one row of each detection group that stands for it.

    Several rows at one detection are one measurement under several
    orientations, so a track states it once.  The row kept is the one on the
    LARGEST cluster -- the track with the most evidence behind it -- with the
    cluster id and then the row position breaking ties, so the choice is a
    function of the state alone.  ``rows_sel`` restricts the reading to a
    subset of rows, which is what naming a culprit needs."""
    gid = np.asarray(gid, np.int64)
    slot_c = np.asarray(slot_c, np.int64)
    n_obs_pt = np.asarray(n_obs_pt, np.int64)
    pos = (
        np.arange(len(gid), dtype=np.int64)
        if rows_sel is None
        else np.asarray(rows_sel, np.int64)
    )
    g = gid[pos]
    c = slot_c[pos]
    ordered = pos[np.lexsort((pos, c, -n_obs_pt[c], g))]
    go = gid[ordered]
    first = np.concatenate(([True], go[1:] != go[:-1]))
    return go[first], ordered[first]


def solve_tracks(cam, state, uv, slot_i, tracks, floor_rad):
    """``(xyzw, verdicts, median reprojection)`` over explicit row sets.

    One :func:`seed_relax.structure.estimate_points_verdicts` call over every
    track, each track's rows owning one point, floored at the member's own
    angular tolerance with the cheirality refusal on."""
    if not tracks:
        return np.zeros((0, 4)), np.zeros(0, np.uint8), np.zeros(0)
    flat = np.concatenate(tracks)
    sizes = np.array([len(t) for t in tracks], np.int64)
    owner = np.repeat(np.arange(len(tracks), dtype=np.int64), sizes)
    pts, at_inf, _census, verdicts, _pruned = structure.estimate_points_verdicts(
        cam,
        state["quats"],
        state["trans"],
        uv[flat],
        slot_i[flat],
        owner,
        len(tracks),
        floor_rad,
    )
    rot, cen = structure.centres_of(state)
    res = structure.reprojection(
        cam, rot, cen, pts, at_inf, uv[flat], slot_i[flat], owner
    )
    med = np.full(len(tracks), np.nan)
    ends = np.cumsum(sizes)
    for k, (lo, hi) in enumerate(zip(np.concatenate(([0], ends[:-1])), ends)):
        v = res[lo:hi]
        v = v[np.isfinite(v)]
        if len(v):
            med[k] = float(np.median(v))
    xyzw = np.concatenate([pts, np.where(at_inf, 0.0, 1.0)[:, None]], axis=1)
    return xyzw, verdicts, med


def widest_angles(cam, state, uv, slot_i, slot_c, want, n_points):
    """Each wanted point's widest ray pair in radians; NaN where it has none."""
    rot, _cen = structure.centres_of(state)
    d = structure.world_rays(cam, rot, uv, slot_i)
    order = np.argsort(slot_c, kind="stable")
    bounds = np.searchsorted(
        np.asarray(slot_c, np.int64)[order], np.arange(int(n_points) + 1)
    )
    theta = np.full(int(n_points), np.nan)
    for p in np.asarray(want, np.int64):
        dk = d[order[bounds[p] : bounds[p + 1]]]
        if len(dk) >= 2:
            theta[p] = math.acos(max(-1.0, min(1.0, float((dk @ dk.T).min()))))
    return theta


def knot_readings(
    cam, state, uv, slot_i, slot_c, gid, shared, knot, n_knots, tol_rad, tol_px
):
    """What each knot's union track says, against what the state already held.

    Per knot: the knot's own rows, its union track (one row per detection), the
    union's verdict and median reprojection, the same rows' median under the
    separate points floored at the member's own tolerance in pixels, the
    largest member distance to the union in units of that member's own depth
    uncertainty, and whether the union both solves and passes.

    ``ok`` is the merge test: the union solves finite, it explains the union's
    own rows at least as well as the separate points did, and every finite
    member sits within ``depth * tol_rad / sin(theta)`` of it.  ``thin`` is the
    union that measures a direction and no depth."""
    n_pts = len(np.asarray(state["clusters"]))
    n_obs_pt = np.bincount(slot_c, minlength=n_pts)
    _g, rep = representatives(gid, slot_c, n_obs_pt)
    rep_of_gid = np.zeros(int(gid.max()) + 1 if len(gid) else 1, np.int64)
    rep_of_gid[gid[rep]] = rep

    kr = np.nonzero(knot[slot_c] >= 0)[0]
    kr = kr[np.argsort(knot[slot_c[kr]], kind="stable")]
    kb = np.searchsorted(knot[slot_c[kr]], np.arange(n_knots + 1))
    knot_rows = [kr[kb[k] : kb[k + 1]] for k in range(n_knots)]
    knot_points = [np.unique(slot_c[r]) for r in knot_rows]
    union_rows = [np.unique(rep_of_gid[gid[r]]) for r in knot_rows]

    xyzw, verdicts, union_med = solve_tracks(
        cam, state, uv, slot_i, union_rows, tol_rad
    )
    rot, cen = structure.centres_of(state)
    own = structure.reprojection(
        cam,
        rot,
        cen,
        np.asarray(state["points"], float),
        np.asarray(state["at_inf"], bool),
        uv,
        slot_i,
        slot_c,
    )
    theta = widest_angles(
        cam, state, uv, slot_i, slot_c, np.nonzero(knot >= 0)[0], n_pts
    )
    pts_state = np.asarray(state["points"], float)
    inf_state = np.asarray(state["at_inf"], bool)

    ok = np.zeros(n_knots, bool)
    thin = np.zeros(n_knots, bool)
    ref_px = np.full(n_knots, np.nan)
    own_med = np.full(n_knots, np.nan)
    z_max = np.full(n_knots, np.nan)
    for k in range(n_knots):
        u = union_rows[k]
        o = own[u]
        o = o[np.isfinite(o)]
        own_med[k] = float(np.median(o)) if len(o) else np.nan
        ref_px[k] = max(own_med[k] if len(o) else 0.0, tol_px)
        thin[k] = int(verdicts[k]) == THIN_VERDICT
        if xyzw[k, 3] == 1.0:
            sh = u[shared[u]]
            c0 = cen[slot_i[int(sh[0]) if len(sh) else int(u[0])]]
            vals = []
            for p in knot_points[k]:
                if inf_state[p] or not np.isfinite(theta[p]) or theta[p] <= 0.0:
                    continue
                dep = float(np.linalg.norm(pts_state[p] - c0))
                unc = dep * float(tol_rad) / max(math.sin(theta[p]), 1e-9)
                vals.append(
                    float(np.linalg.norm(pts_state[p] - xyzw[k, :3]) / max(unc, 1e-12))
                )
            if vals:
                z_max[k] = float(max(vals))
        ok[k] = bool(
            int(verdicts[k]) in FINITE_VERDICTS
            and np.isfinite(union_med[k])
            and union_med[k] <= ref_px[k]
            and (not np.isfinite(z_max[k]) or z_max[k] <= 1.0)
        )
    return {
        "n_obs_pt": n_obs_pt,
        "knot_rows": knot_rows,
        "knot_points": knot_points,
        "union_rows": union_rows,
        "xyzw": xyzw,
        "verdicts": verdicts,
        "union_med": union_med,
        "own_med": own_med,
        "ref_px": ref_px,
        "z_max": z_max,
        "ok": ok,
        "thin": thin,
    }


def drop_one(cam, state, uv, slot_i, slot_c, shared, rd, tol_rad):
    """``{knot: (culprit, reprojection, solved point)}`` where one removal helps.

    For every knot whose union contradicts, each member's EXCLUSIVE detections
    are dropped in turn and the rest resolved.  Where a removal restores a
    finite union under the same bar the union was judged on, that member is
    what the arithmetic votes against; the best such removal is taken, with the
    point id breaking a tie.  This is drop-one-and-resolve, the shape the
    cheirality-minority prune already uses."""
    tracks, owners = [], []
    for k in np.nonzero(~rd["ok"] & ~rd["thin"])[0]:
        u = rd["union_rows"][k]
        su = shared[u]
        for p in rd["knot_points"][k]:
            drop = (slot_c[u] == p) & ~su
            keep = u[~drop]
            if len(keep) >= 2 and drop.any():
                tracks.append(keep)
                owners.append((int(k), int(p)))
    xyzw, verdicts, med = solve_tracks(cam, state, uv, slot_i, tracks, tol_rad)
    best = {}
    for t, (k, p) in enumerate(owners):
        if int(verdicts[t]) not in FINITE_VERDICTS or not np.isfinite(med[t]):
            continue
        if med[t] > rd["ref_px"][k]:
            continue
        cur = best.get(k)
        if cur is None or med[t] < cur[1] or (med[t] == cur[1] and p < cur[0]):
            best[k] = (p, float(med[t]), xyzw[t])
    return best


# -------------------------------------------------------------------- the stage


def reconcile_points(mx, cam, state, refine_radius, tol_rad, f_eq, trace=None):
    """``(member, state, census)`` with every knot the arithmetic settles.

    The member that comes back carries the merged tracks: an absorbed point's
    rows name the surviving cluster, the duplicate rows on one detection leave
    the admission, and a culled member's rows leave it too.  A knot no single
    removal reconciles is left exactly as it stood, and counted.

    The stage REFUSES, and changes nothing, where the state holds no
    observations, where the member states no affine shapes, or where the source
    file states no refine radius: the unit scale the relation is stated in is
    then unstated, and a guess at it would merge tracks on a number nobody
    measured."""
    census = {
        "unit_frac": float(NEAR_UNIT_FRAC),
        "radius_frac": float(NEAR_RADIUS_FRAC),
    }
    rows, slot_i, slot_c = structure.state_rows(mx, state)
    clusters = np.asarray(state["clusters"], np.int64)
    census["n_rows"] = int(len(rows))
    census["n_points"] = int(len(clusters))
    if not len(rows):
        census["refused"] = "the state holds no observations"
        return mx, state, census
    if mx.obs_shape is None:
        census["refused"] = "the member states no affine shapes"
        return mx, state, census
    if not refine_radius:
        census["refused"] = "the source file states no refine radius"
        return mx, state, census

    uv = np.ascontiguousarray(np.asarray(mx.obs_uv, float)[rows])
    r_px = evict.feature_radius(np.asarray(mx.obs_shape, float)[rows], refine_radius)
    gid = detection_groups(slot_i, uv)
    pair_a, pair_b = near_rows(slot_i, slot_c, uv, r_px, refine_radius)
    knot, n_knots, shared = knot_components(slot_c, gid, pair_a, pair_b, len(clusters))
    census.update(
        {
            "n_detections": int(gid.max() + 1),
            "n_near_pairs": int(len(pair_a)),
            "n_obs_shared": int(shared.sum()),
            "n_points_in_knot": int((knot >= 0).sum()),
            "n_knots": int(n_knots),
            "merged": 0,
            "merged_bearing": 0,
            "culled": 0,
            "refused_knots": 0,
            "rows_deduped": 0,
            "n_points_dropped": 0,
            "n_points_after": int(len(clusters)),
        }
    )
    if not n_knots:
        if trace is not None:
            trace("    reconcile: no knot, nothing rests on one measurement")
        return mx, state, census

    tol_px = float(tol_rad) * float(f_eq)
    census["tol_px"] = tol_px
    rd = knot_readings(
        cam, state, uv, slot_i, slot_c, gid, shared, knot, n_knots, tol_rad, tol_px
    )
    best = drop_one(cam, state, uv, slot_i, slot_c, shared, rd, tol_rad)

    new_obs_c = np.asarray(mx.obs_c, np.int64).copy()
    drop_pt = np.zeros(len(clusters), bool)
    prune = np.zeros(len(rows), bool)
    points = np.asarray(state["points"], float).copy()
    at_inf = np.asarray(state["at_inf"], bool).copy()
    n_obs_pt = rd["n_obs_pt"]
    n_merge = n_bearing = n_cull = n_refuse = 0
    for k in range(n_knots):
        culprit = None
        if rd["thin"][k] or rd["ok"][k]:
            survivors = rd["knot_points"][k]
            keep_rows = rd["union_rows"][k]
            solved = rd["xyzw"][k]
        elif k in best:
            culprit = best[k][0]
            survivors = rd["knot_points"][k][rd["knot_points"][k] != culprit]
            sel = rd["knot_rows"][k][slot_c[rd["knot_rows"][k]] != culprit]
            _g, reps = representatives(gid, slot_c, n_obs_pt, sel)
            keep_rows = np.sort(reps)
            solved = best[k][2]
        else:
            n_refuse += 1
            continue
        if not len(survivors) or not len(keep_rows):
            n_refuse += 1
            continue
        head = int(survivors[np.lexsort((survivors, -n_obs_pt[survivors]))[0]])
        mine = rd["knot_rows"][k]
        kept = np.isin(mine, keep_rows)
        new_obs_c[rows[keep_rows]] = int(clusters[head])
        prune[mine[~kept]] = True
        drop_pt[survivors] = True
        drop_pt[head] = False
        points[head] = solved[:3]
        at_inf[head] = solved[3] == 0.0
        if rd["thin"][k]:
            n_bearing += 1
        elif rd["ok"][k]:
            n_merge += 1
        else:
            n_cull += 1
            drop_pt[culprit] = True

    census.update(
        {
            "merged": int(n_merge),
            "merged_bearing": int(n_bearing),
            "culled": int(n_cull),
            "refused_knots": int(n_refuse),
            "rows_deduped": int(prune.sum()),
            "n_points_dropped": int(drop_pt.sum()),
            "n_points_after": int(len(clusters) - drop_pt.sum()),
        }
    )
    if trace is not None:
        trace(
            f"    reconcile: {n_knots} knots over {census['n_points_in_knot']} "
            f"points, {n_merge} merged, {n_bearing} as bearings, {n_cull} "
            f"culled, {n_refuse} refused, {census['rows_deduped']} rows deduped"
        )
    if not drop_pt.any() and not prune.any():
        return mx, state, census

    ev = _ev()
    mx = ev.Member(
        mx.idx,
        mx.model,
        list(mx.names),
        mx.camera,
        mx.f_eq,
        mx.rvec,
        mx.tvec,
        mx.posed,
        mx.pts,
        (new_obs_c, mx.obs_i, mx.obs_uv, mx.obs_f),
        shapes=mx.obs_shape,
        keep=mx.keep_mask(),
        dropped=mx.dropped,
    )
    state = structure.with_pruned(state, rows[prune])
    keep_pt = ~drop_pt
    state["clusters"] = clusters[keep_pt]
    state["points"] = points[keep_pt]
    state["at_inf"] = at_inf[keep_pt]
    return mx, state, census


def reconcile_stage(mx, cam, state, refine_radius, tol_rad, f_eq, trace=None):
    """The stage as the pipeline runs it: relation, readings, resolve, apply."""
    mx, state, census = reconcile_points(
        mx, cam, state, refine_radius, tol_rad, f_eq, trace
    )
    if trace is not None and census.get("refused"):
        trace(f"    reconcile: refused ({census['refused']})")
    return mx, state, census
