# Spike: mutual-kNN edges/clusters as post-reconstruction track candidates

Prototype for the design direction in `specs/core/mutual-knn-matching.md`
("Opportunity: mutual edges as post-reconstruction track candidates"): validate
mutual-kNN candidates against a *known* cluster-based reconstruction (cheap
geometry) instead of blind two-view RANSAC. **Throwaway measurement code.**

Two harnesses, and the second corrects the first:

- `spike.py` — tests mutual **edges** individually (two-view triangulation).
- `cluster_spike.py` — tests mutual **clusters** (connected components of the
  edges = putative multi-view tracks), N-view triangulated as a unit. **This is
  the honest unit; trust it over the edge numbers.**

Run: `pixi run -e test python prototypes/mutual_edge_densify/{spike,cluster_spike}.py`
(expects cluster reconstructions staged at `/tmp/<dataset>_ab_ws/cluster_recon.sfmr`).

## Headline

Two conclusions survive; one early number does not.

1. **Cost thesis holds (strongly).** Validating candidates against known poses
   is cheap and RANSAC-free: edge/cluster build ~1–3s, the geometric test
   0.03–0.22s, vs the 16–26s the *standalone* `--mutual-knn` matcher spends on
   RANSAC. Known geometry replaces per-pair model *estimation* with one
   triangulation + reprojection.
2. **Non-overlapping pairs gate for free** (edge-level): overlapping pairs
   accept 14–18%, non-overlapping 1.9–2.4% — the pose-based reprojection test
   rejects the junk with no covisibility precompute.
3. **The yield is much smaller than edge-level suggested.** The edge spike
   implied ~12k "new tracks" on seattle. At cluster granularity that collapses
   to a few **hundred** candidate tracks / a few percent more observations —
   because the raw mutual edges form one giant over-merged component, and most
   "accepted edges" lived *inside* it without forming coherent tracks.

## Why edge-level overcounts: catastrophic over-merge

Taken as a graph, the raw mutual edges (triangle filter off) connect almost
every feature into a **single** component:

| dataset | edges | nodes | clusters | largest cluster |
|---|---|---|---|---|
| seattle_backyard | 363,575 | 83,033 | **41** | **82,946** (≈ all nodes) |
| seoul_bull | 167,750 | 36,859 | **14** | 36,828 (≈ all nodes) |

This is the connected-components pathology the spec warns about. The geometric
test correctly *rejects* the blob (0 of the ≥10-member clusters accepted at
tmin=0), so edge-level "accepted edges" (56k on seattle) are mostly intra-blob
edges that pass pairwise reprojection but do **not** form valid tracks. Counting
edges overstates yield by 1–2 orders of magnitude.

## Cluster-level yield vs the triangle filter (`cluster_spike.py`)

The triangle filter (require an edge be corroborated by `triangle_min` shared
mutual neighbours) shatters the blob. seattle_backyard:

| triangle_min | clusters | largest | accepted tracks | new | extend | merge | ~new obs |
|---|---|---|---|---|---|---|---|
| 0 | 41 | 82,946 | 9 | 8 | 1 | 0 | ~16 |
| 1 | 435 | 65,075 | 67 | 13 | 53 | 1 | ~49 |
| 2 | 1,701 | 44,191 | 434 | 193 | 238 | 3 | ~444 |
| 4 | 3,018 | 9,563 | 1,140 | 295 | 825 | 20 | ~752 |

seoul_bull tops out around ~270 new observations (tmin=4). Pose self-check was
sub-pixel (0.23 / 0.28px) throughout, so the geometry is sound.

Read against the ~20k existing observations in the seattle reconstruction, the
best case is a **few-percent** densification — real, but a long way from the
edge-level impression.

## Honest caveats

- **Over-merge is the central unsolved problem, not a detail.** Even at tmin=4
  a multi-thousand-member blob survives, `clusters≥10 accepted` and the
  `dup-image` count (clusters with two features in the *same* image — a
  track-consistency violation) both grow. Connected-components + triangle filter
  is not enough; a real feature needs per-cluster splitting / dedup (one feature
  per image), not just CC.
- **Yield ≠ reconstruction improvement.** These are *candidate* tracks. Whether
  adding them improves density **and accuracy** needs an actual add +
  re-triangulate + bundle-adjust + measure — not done here.
- **Geometric consistency is necessary, not sufficient.** Reprojection agreement
  with known poses can still admit wrong matches on repeated structure.
- **`triangle_min` and the 4px / 1.5° gates are unswept knobs**; the right
  `triangle_min` is dataset-dependent (higher = fewer, cleaner clusters).

## Verdict

The reframe is mechanically sound and *cheap*, and your instinct to evaluate
clusters rather than edges was the right correction — it cut an inflated
~12k-track number down to the few-hundred reality. But as it stands the
opportunity is **modest and gated by an over-merge problem** that connected
components alone doesn't solve. Worth pursuing only with (a) a real cluster
splitter and (b) an end-to-end density+accuracy measurement after BA. Otherwise
the standalone matcher and this densification idea are both narrow wins at best.

## If it graduates

1. Emit mutual edges as a byproduct of `--cluster` (shared k-NN query).
2. Cluster with the triangle filter **plus** a splitter enforcing ≤1 feature per
   image and geometric coherence per component.
3. N-view triangulate, trim, reconcile with existing tracks, add, re-triangulate,
   bundle-adjust.
4. Measure final density **and accuracy** on seattle (win), dino (standalone
   loss), kerry, seoul.
