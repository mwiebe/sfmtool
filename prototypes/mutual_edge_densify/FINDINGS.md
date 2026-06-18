# Spike: mutual-kNN edges as post-reconstruction track candidates

Prototype validating the design direction recorded in
`specs/core/mutual-knn-matching.md` ("Opportunity: mutual edges as
post-reconstruction track candidates"). **Throwaway measurement code**, not a
feature. Run: `pixi run -e test python prototypes/mutual_edge_densify/spike.py`
(expects cluster reconstructions staged under `/tmp/<dataset>_ab_ws/`).

## What it does

1. Loads a **cluster-based** reconstruction (the cheap, high-precision graph).
2. **Self-checks the pose math** by reprojecting existing tracks — must be
   sub-pixel or the run aborts.
3. Computes the **raw mutual edges** (`mutual_knn_matches`, pre-verification).
4. Tests every edge against *known* geometry: triangulate from the two known
   camera centers, require cheirality (in front of both), a real triangulation
   angle (≥1.5°), and per-view reprojection error < 4px. **No RANSAC.**
5. Reports pass rate, accept rate on overlapping vs non-overlapping pairs, the
   densification yield, and cost.

## Results (2026-06)

Pose self-check passed both datasets (median reproj 0.23px / 0.28px).

| dataset | candidates | edge compute | geom test (no RANSAC) | standalone RANSAC¹ |
|---|---|---|---|---|
| seattle_backyard | 363,575 | 2.87s | **2.77s** | 16.5s |
| seoul_bull | 167,750 | 0.97s | **0.84s** | 26.4s |

¹ the standalone `--mutual-knn` matcher's verification time on the same edges.

**Free gating of non-overlapping pairs** (the standalone matcher's cost sink):

| dataset | overlapping pairs accept | non-overlapping pairs accept |
|---|---|---|
| seattle_backyard | 18.5% (296,669 cand) | **1.9%** (66,906 cand) |
| seoul_bull | 14.2% (72,297 cand) | **2.4%** (95,453 cand) |

**Densification yield** (accepted edges, by track membership of endpoints):

| dataset | new-track | extend-track | merge | already-linked |
|---|---|---|---|---|
| seattle_backyard | 12,136 | 10,728 | 4,454 | 28,706 |
| seoul_bull | 8,314 | 1,428 | 123 | 2,639 |

(seattle reconstruction has ~20k existing observations, so the new+extend edges
are a large fraction of the existing structure.)

## Verdict

The hypothesis holds on both datasets:

- **Cost flips.** Validating the full raw edge set against known poses is ~6×
  cheaper than the RANSAC the standalone matcher spends, because known geometry
  replaces per-pair model *estimation* with one triangulation + two
  reprojections.
- **Non-overlapping pairs gate for free.** 1.9–2.4% accept vs 14–18% on
  overlapping pairs — the pose-based reprojection test rejects the junk with no
  covisibility precompute.
- **There is real densification headroom** the cluster reconstruction misses.

## Honest caveats (why this is a spike, not a result to ship on)

- **"new-track"/"extend" counts are edges, not tracks** — an upper bound on
  yield. Real yield needs connected-component track assembly from accepted
  edges (multiple edges collapse into one track) plus re-triangulation.
- **Geometric consistency ≠ correct match.** Reprojection agreement with known
  poses is necessary, not sufficient: repeated structure can produce a
  consistent-but-wrong 3D point. Triangulation+cheirality+angle is stronger than
  a pure epipolar test, but a real implementation must gate adds with bundle
  adjustment and outlier rejection before trusting accuracy.
- **Threshold sensitivity.** Yield depends on the 4px / 1.5° gates; not swept
  here. The overlap-vs-non-overlap *ratio* is the robust finding.
- **The Python edge-gather loop dominates the geom-test time;** a real
  implementation in the Rust core would be far faster still.

## If this graduates to a feature

1. Emit raw mutual edges as a near-free byproduct of the `--cluster` run (shared
   k-NN query).
2. After the initial solve, run the geometric gate in the Rust core.
3. Assemble accepted edges into connected-component tracks, add observations /
   seed points, re-triangulate, and bundle-adjust.
4. Measure final reconstruction **density and accuracy** (not just edge counts)
   on seattle (win), dino (standalone-loss), kerry, seoul.
