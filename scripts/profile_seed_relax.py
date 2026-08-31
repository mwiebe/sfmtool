# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Where the seed relaxation spends its time, on one stored member.

The relaxation runs eight stages on a member and its source selection handle
(`scripts/seed_relax/pipeline.py`): gate, relax, fill in, reconcile, hand over,
release, re-estimate, report.  This harness runs that chain once on a member
restored from a study sidecar and reports a wall-clock table: the eight stages,
and inside each the pieces it is made of -- the fill-in's join, ring assignment,
point estimation and held adjustment; the reconciliation's tangle building,
union solves and drop-one resolve; the hand-over's containment rule,
two-observation sweep and its own held adjustment.  Both the reconciliation and
the hand-over read their candidate pairs from the same keypoint-reach kernel,
which is timed under each of them separately.

The timing is taken by wrapping the package's own functions where the caller
looks them up, so the chain runs exactly as it ships and nothing in the package
reads a clock.  A name several stages call -- the held adjustment, the keypoint
reach kernel -- is attributed to the stage that was running, so one row is one
stage's use of it.

Usage::

    python scripts/profile_seed_relax.py CANDIDATE_SOLVES_DIR IDX MATCHES_PATH

``CANDIDATE_SOLVES_DIR`` holds the ``member_arrays.npz`` sidecar and its
``manifest.json``; ``IDX`` is the hypothesis index of a rotation-only member in
that sidecar; ``MATCHES_PATH`` is the capture's cluster-patches ``.matches``
file, opened the way the run opens it.  It must be the file the member was
drawn from: the fill-in joins member rows to file rows by ``(image, feature)``
and refuses a handle whose image table is not the member's, so a
seed-restricted file -- which names only the seed's own images -- leaves the
fill-in, and everything downstream of it, unexercised.  ``--repeat N`` runs the
chain N times and reports the median of each row.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import zlib
from collections import defaultdict
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parent
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

#: The minimum cluster span the loader admits, matching `exp_fast_seed`.
MIN_SPAN = 2

#: Separates a stage from a name it called, in a scoped clock key.
SCOPE = " :: "


def load_member(cs_dir, idx):
    """One rotation-only member of a study sidecar, as it was committed."""
    import seed_candidate_eval as EV

    d = np.load(Path(cs_dir) / "member_arrays.npz")
    meta = json.loads(zlib.decompress(d["_meta"]).decode("utf-8"))
    key = f"m{int(idx):04d}"
    arr = {k.split("__", 1)[1]: d[k] for k in d.files if k.startswith(key + "__")}
    mm = meta["members"][key]
    return EV.member_from_arrays(
        {
            "idx": int(idx),
            "model": mm["model"],
            "names": meta["names"],
            "camera": mm["camera"],
            "f_eq": mm["f_eq"],
            "rvec": arr["rvec"],
            "tvec": arr["tvec"],
            "posed": arr["posed"],
            "pts": arr.get("pts"),
            "obs_c": arr["obs_c"],
            "obs_i": arr["obs_i"],
            "obs_uv": arr["obs_uv"],
            "obs_f": arr.get("obs_f"),
            "obs_shape": arr.get("obs_shape"),
            "keep": arr.get("keep"),
        }
    )


def open_source(matches_path):
    """The selection handle the run holds, from the capture's matches file."""
    from sfmtool._sfmtool.io import MatchesFile

    return MatchesFile(str(matches_path)).select_clusters(min_span=MIN_SPAN)


def source_mismatch(m, source):
    """Why this handle cannot fill this member in, or ``None``.

    The fill-in's join is by ``(image, feature)`` against the file's own image
    table, so it refuses a handle whose names are not the member's -- which is
    what a seed-restricted file is, and what leaves every stage from the
    fill-in onward doing nothing."""
    names = [str(n).replace("\\", "/") for n in source.image_names]
    if names == list(m.names):
        return None
    return (
        f"the matches file names {len(names)} images and the member {len(m.names)}: "
        "the fill-in reads its clusters from the capture's own cluster-patches "
        "file and refuses a handle whose image table is not the member's, so "
        "the fill-in and every stage after it will do nothing"
    )


class Clock:
    """Wall-clock totals and call counts, by name and by calling stage."""

    def __init__(self):
        self.total = defaultdict(float)
        self.calls = defaultdict(int)
        self.each = defaultdict(list)
        self.stage = None

    def _record(self, label, dt):
        self.total[label] += dt
        self.calls[label] += 1
        self.each[label].append(dt)

    def wrap(self, module, name, label=None, stage=False):
        """Time every call of ``module.name`` where its callers look it up.

        ``stage=True`` marks the call a stage boundary: everything timed
        underneath it is recorded a second time under ``<stage> :: <name>``, so
        a name two stages call reads as one row per stage."""
        if label is None:
            label = f"{module.__name__.rsplit('.', 1)[-1]}.{name}"
        fn = getattr(module, name)

        def timed(*a, **kw):
            outer = self.stage
            if stage:
                self.stage = label
            t0 = time.perf_counter()
            try:
                return fn(*a, **kw)
            finally:
                dt = time.perf_counter() - t0
                self.stage = outer
                self._record(label, dt)
                if not stage and outer is not None:
                    self._record(outer + SCOPE + label, dt)

        timed.__name__ = getattr(fn, "__name__", name)
        timed.__doc__ = getattr(fn, "__doc__", None)
        setattr(module, name, timed)
        return fn

    def unscoped(self, label):
        """``(seconds, calls)`` of the calls no stage boundary was inside."""
        secs = self.total.get(label, 0.0)
        calls = self.calls.get(label, 0)
        for key, v in self.total.items():
            if key.endswith(SCOPE + label):
                secs -= v
                calls -= self.calls[key]
        return secs, calls


#: The stage boundaries, in the order the chain runs them.
STAGES = [
    ("lens", "rot_lens_ba", "stage 1 lens on bearings"),
    ("relaxation", "relax_oriented", "stage 2 relaxation"),
    ("fill", "fill_in", "stage 3 fill-in"),
    ("reconcile", "reconcile_stage", "stage 4 reconcile"),
    ("evict", "evict_stage", "stage 5 hand-over"),
    ("lens", "release_at_knots", "stage 6 late lens release"),
    ("report", "runaway_report", "stage 8 runaway report"),
]

#: The keypoint-reach kernel, timed where `seed_relax.pairs` looks it up.
REACH_KERNEL = "pairs.keypoint_pairs_within_reach"

#: Rows of the table: the label, and how deep it sits under its parent.
ROWS = [
    ("stage 1 lens on bearings", 0),
    ("stage 2 relaxation", 0),
    ("graph.member_graph", 1),
    ("graph.stage_pairs", 1),
    ("graph.pair_rays", 2),
    ("graph.baseline_direction", 2),
    ("averaging.centres_by_averaging", 1),
    ("orientation.angw_bit", 1),
    ("structure.triangulate_placed", 1),
    ("structure.grow_more", 1),
    ("structure.build_ba_inputs", 1),
    ("structure.stage_adjust", 1),
    ("stage 3 fill-in", 0),
    ("fill.source_clusters", 1),
    ("rings.band_order", 1),
    ("fill.extend_member", 1),
    ("fill.ring_rows", 1),
    ("fill.estimate_points", 1),
    ("stage 3 fill-in :: fill.adjust_held", 1),
    ("stage 4 reconcile", 0),
    ("reconcile.detection_groups", 1),
    ("reconcile.near_rows", 1),
    ("stage 4 reconcile :: " + REACH_KERNEL, 2),
    ("reconcile.tangle_components", 1),
    ("reconcile.tangle_readings", 1),
    ("reconcile.representatives", 2),
    ("reconcile.widest_angles", 2),
    ("reconcile.solve_tracks", 2),
    ("reconcile.drop_one", 1),
    ("stage 5 hand-over", 0),
    ("evict.covered_by_finer", 1),
    ("stage 5 hand-over :: " + REACH_KERNEL, 2),
    ("evict.two_observation_sweep", 1),
    ("evict.band_census", 1),
    ("rings.assign_rings", 2),
    ("stage 5 hand-over :: fill.adjust_held", 1),
    ("stage 6 late lens release", 0),
    ("stage 7 re-estimate", 0),
    ("structure.state_rows (final)", 1),
    ("structure.estimate_points (final)", 1),
    ("structure.reprojection (final)", 1),
    ("stage 8 runaway report", 0),
    # A ninth stage is wrapped here as soon as the chain grows one.
    ("run_member", 0),
]


def instrument(clock):
    """Wrap every timed name, in the module its caller reads it from."""
    import seed_relax
    from seed_relax import (
        averaging,
        evict,
        fill,
        graph,
        lens,
        orientation,
        pipeline,
        reconcile,
        relaxation,
        report,
        rings,
        structure,
    )

    mods = {
        "averaging": averaging,
        "evict": evict,
        "fill": fill,
        "graph": graph,
        "lens": lens,
        "orientation": orientation,
        "reconcile": reconcile,
        "relaxation": relaxation,
        "report": report,
        "rings": rings,
        "structure": structure,
    }
    for mod, name, label in STAGES:
        clock.wrap(mods[mod], name, label, stage=True)

    clock.wrap(graph, "member_graph")
    clock.wrap(graph, "stage_pairs")
    clock.wrap(graph, "pair_rays")
    clock.wrap(graph, "baseline_direction")
    clock.wrap(averaging, "centres_by_averaging")
    clock.wrap(orientation, "angw_bit")
    clock.wrap(structure, "triangulate_placed")
    clock.wrap(structure, "grow_more")
    clock.wrap(structure, "build_ba_inputs")
    clock.wrap(structure, "stage_adjust")
    clock.wrap(structure, "estimate_points")
    clock.wrap(structure, "reprojection")
    clock.wrap(structure, "state_rows")
    clock.wrap(fill, "source_clusters")
    clock.wrap(fill, "extend_member")
    clock.wrap(fill, "ring_rows")
    clock.wrap(fill, "adjust_held")
    # `fill` binds the estimator by name at import, so the fill-in's own calls
    # are counted separately from the re-estimation's.
    clock.wrap(fill, "estimate_points")
    clock.wrap(rings, "assign_rings")
    clock.wrap(rings, "band_order")
    clock.wrap(reconcile, "detection_groups")
    clock.wrap(reconcile, "near_rows")
    clock.wrap(reconcile, "tangle_components")
    clock.wrap(reconcile, "tangle_readings")
    clock.wrap(reconcile, "representatives")
    clock.wrap(reconcile, "widest_angles")
    clock.wrap(reconcile, "solve_tracks")
    clock.wrap(reconcile, "drop_one")
    clock.wrap(evict, "covered_by_finer")
    clock.wrap(evict, "two_observation_sweep")
    clock.wrap(evict, "band_census")
    # `seed_relax.pairs` imports the kernel inside the call, so the wrap has to
    # sit on the binding module itself rather than on a name `pairs` holds.
    import sfmtool._sfmtool.analysis as analysis

    clock.wrap(analysis, "keypoint_pairs_within_reach", REACH_KERNEL)
    return seed_relax, pipeline


def run_once(cs_dir, idx, matches_path, source):
    """One whole chain, timed.  Returns the clock and the result."""
    clock = Clock()
    seed_relax, pipeline = instrument(clock)
    m = load_member(cs_dir, idx)
    opts = seed_relax.Options()
    t0 = time.perf_counter()
    result = pipeline.run_member(m, source, opts)
    clock.total["run_member"] = time.perf_counter() - t0
    clock.calls["run_member"] = 1
    # The re-estimation is the pipeline's own, written inline: what is left of
    # these three names once the calls made inside a stage are taken out.
    for name in ("state_rows", "estimate_points", "reprojection"):
        secs, calls = clock.unscoped(f"structure.{name}")
        clock.total[f"structure.{name} (final)"] = secs
        clock.calls[f"structure.{name} (final)"] = calls
    clock.total["stage 7 re-estimate"] = sum(
        clock.total[f"structure.{n} (final)"]
        for n in ("state_rows", "estimate_points", "reprojection")
    )
    clock.calls["stage 7 re-estimate"] = 1
    return clock, result


def render(clock, total):
    """The table, one row per timed name."""
    lines = [f"{'stage':<52}{'calls':>8}{'seconds':>12}{'% of run':>10}"]
    lines.append("-" * 82)
    for label, depth in ROWS:
        if label not in clock.total:
            continue
        secs = clock.total[label]
        name = ("  " * depth) + label.split(SCOPE)[-1]
        pct = 100.0 * secs / total if total > 0 else 0.0
        lines.append(f"{name:<52}{clock.calls[label]:>8}{secs:>12.3f}{pct:>10.1f}")
    return "\n".join(lines)


def per_ring(clock):
    """The fill-in's per-ring rows, in the order the rings ran."""
    names = [
        ("fill.ring_rows", "ring_rows"),
        ("fill.estimate_points", "estimate_points"),
        ("stage 3 fill-in :: fill.adjust_held", "adjust_held"),
    ]
    n = max((len(clock.each[k]) for k, _ in names), default=0)
    if not n:
        return "no ring ran"
    lines = [f"{'ring':<8}" + "".join(f"{short:>22}" for _k, short in names)]
    for r in range(n):
        row = f"{r:<8}"
        for k, _short in names:
            v = clock.each[k]
            row += f"{(v[r] if r < len(v) else float('nan')):>22.3f}"
        lines.append(row)
    for k in ("fill.source_clusters", "rings.band_order", "fill.extend_member"):
        lines.append(f"{k:<30}{clock.total[k]:>10.3f} s")
    return "\n".join(lines)


def say_census(result):
    """What each stage after the relaxation decided, or why it refused."""
    fc = result.census.get("fill", {})
    if fc.get("refused"):
        print(f"  fill-in: refused ({fc['refused']})")
    else:
        print(
            f"  fill-in: {fc.get('n_candidates')} candidates in "
            f"{fc.get('n_rings')} rings, {fc.get('n_added')} added"
        )
    rc = result.census.get("reconcile", {})
    if rc.get("refused") or rc.get("held"):
        print(f"  reconcile: {rc.get('refused') or 'held by ' + rc['held']}")
    else:
        print(
            f"  reconcile: {rc.get('n_tangles')} tangles, {rc.get('merged')} merged, "
            f"{rc.get('culled')} culled, {rc.get('refused_tangles')} refused"
        )
    ec = result.census.get("evict", {})
    if ec.get("refused") or ec.get("held"):
        print(f"  hand-over: {ec.get('refused') or 'held by ' + ec['held']}")
    else:
        print(
            f"  hand-over: {ec.get('n_obs_evicted')} observations retired of "
            f"{ec.get('n_rows')}, {ec.get('n_points_kept')} points kept of "
            f"{ec.get('n_points')}"
        )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("candidate_solves", help="directory holding member_arrays.npz")
    ap.add_argument("idx", type=int, help="hypothesis index of the member")
    ap.add_argument("matches", help="the capture's cluster-patches .matches file")
    ap.add_argument("--repeat", type=int, default=1, help="runs to take the median of")
    args = ap.parse_args(argv)

    source = open_source(args.matches)
    why = source_mismatch(load_member(args.candidate_solves, args.idx), source)
    if why is not None:
        print(f"WARNING: {why}", flush=True)

    runs = []
    for r in range(max(1, int(args.repeat))):
        clock, result = run_once(args.candidate_solves, args.idx, args.matches, source)
        runs.append(clock)
        if r == 0:
            state = result.state
            n = 0 if state is None else len(state["at_inf"])
            fin = (
                0 if state is None else int((~np.asarray(state["at_inf"], bool)).sum())
            )
            print(
                f"member h{args.idx}: refused={result.refused} "
                f"points={n} finite={fin} "
                f"placed={0 if state is None else len(state['frames'])}"
            )
            say_census(result)

    merged = Clock()
    for label in {k for c in runs for k in c.total}:
        vals = sorted(c.total.get(label, 0.0) for c in runs)
        merged.total[label] = vals[len(vals) // 2]
        merged.calls[label] = runs[0].calls.get(label, 0)
    print()
    print(render(merged, merged.total.get("run_member", 0.0)))
    print()
    print(per_ring(runs[0]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
