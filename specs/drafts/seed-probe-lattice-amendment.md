# Probe focal on the scan lattice (amendment)

**Status:** Draft

Amends [core/geometry/seed-hypothesis-loop.md](../core/geometry/seed-hypothesis-loop.md)
§ "Choice among seed groups", which specifies the probe focal as the raw
pairwise focal vote and points back here.

## Purpose

Every seed group is probed at one focal before any group is judged. That focal
is the raw vote (under an equidistant context, the verdict's equidistant focal),
so it moves continuously with the vote. The focal scan reads rungs of a lattice
per capture, `max(w, h) * 1.15^k`, so a vote movement under half a step moves no
scanned focal; the probe has no such property. The vote still moves between
files that differ by a fraction of a percent of their members: on `SeoulBull`'s
perturbed cluster files it reads 247.0 to 283.0 px
([measurements](../core/geometry/seed-hypothesis-loop-measurements.md#pick-stability-after-milestone-b-2026-10-09)),
and before the vote read five draws per pair, a probe at 264.7 px instead of
258.3 px grew group 0 to 4 to 7 frames instead of 8 and moved the first
candidate to another seed group
([measurements](../core/geometry/seed-hypothesis-loop-measurements.md#pick-stability-under-small-changes-to-the-cluster-file-2026-10-09)).

## Proposal

The probe runs at the lattice rung nearest the raw focal in log-focal, kept
inside the scan band. With no vote, the nominal `0.9 * max(w, h)` snaps the same
way. Everything that reads the probe focal (the probe, the widen, the outcome
order's median residual) reads the snapped value; the scan window is still
chosen by the raw focal, as it is now.

A vote that moves by less than half a step from its nearest rung then moves
nothing up to the scan. A vote that crosses the log-midpoint between two rungs
still moves the probe by a whole step. Probing at both rungs that bracket the
vote and keeping the better outcome by the outcome order removes that flip, and
doubles the cost of the probe and of every finish that follows it.

Validation is the perturbation suite of the measurements above (13 files per
ground-truth capture) and a fleet run: the snap is adopted when the first
candidate changes on no more files than it does now, no base `h00` regresses
against its ground truth, and the fleet's released focals against their
references do not get worse.

## Open questions

- **Iteration caps.** The bundle adjustments of the probe, the widen and the
  scan stop at fixed evaluation caps (25 to 60). With member order fixed they
  are deterministic on one file, but a capped adjustment still amplifies a
  small change of its input when the result is read against a gate. Whether
  any cap should rise, or the result be read only once converged, is open.
- **Widen-gate margin.** The widen admits a frame when its inlier fraction
  reaches `0.35` times the probe's median per-frame inlier fraction. On
  `KerryPark480` the drop files change which frames the widen admits after the
  rotation core. Whether that gate needs a margin, or a frame near it a second
  reading, is measured by the drop files of the perturbation suite.
