# Mutual-kNN Matching

## The Idea

The background-floor track-cluster matcher (`specs/core/track-cluster-matching.md`)
matches a whole collection by clustering one descriptor corpus instead of
matching image pairs. Its precision mechanism is a **per-descriptor adaptive
radius**: a descriptor keeps cross-image neighbours within `alpha * B_i`, where
`B_i` is its `d`-th-nearest (background) distance.

That radius is excellent for precision but systematically drops the **hard,
wide-baseline matches** that SfM needs most. The cause is structural, not a
tuning miss: a wide-baseline observation of a 3-D point has a genuinely larger
descriptor distance (the viewpoint changed its appearance), so among a
descriptor's nearest neighbours it sits at, typically, **rank ~4** — behind two
or three *short-baseline* views of the same point. The radius `alpha * B_i` is
set by those nearest (short-baseline) neighbours, so the wide one falls just
outside and is cut. The point's track then survives as short-baseline-only:
small parallax, poorly conditioned triangulation, fewer 3-D points.

Mutual-kNN keeps the wide match precisely because it is **rank-bounded, not
radius-bounded**. Each descriptor keeps its `k` nearest *cross-image*
neighbours, and an edge `a-b` survives iff `b` is among `a`'s top-`k` **and**
`a` is among `b`'s — a **mutual** nearest-neighbour. A rank-4 wide match is well
inside the top-`k`, so it is kept; mutuality (not a radius) rejects the random
collisions. The matcher deliberately keeps more, and noisier, candidates and
leans on the downstream **geometric verification** (RANSAC two-view geometry,
already in the pipeline) to reject the false ones — exactly what exhaustive
matching relies on, but driven by one shared approximate k-NN query instead of
all-pairs descriptor matching.

## Empirical Observations

These measurements come from the four bundled datasets (`seoul_bull_sculpture`,
`seattle_backyard`, `kerry_park`, `dino_dog_toy`). "Hard" pairs are the
bottom-tercile image pairs by verified inlier count (the weak / wide-baseline
overlaps); recall is against exhaustive matching's geometrically-verified
inliers.

### The floor's hard-pair recall deficit is systematic

| dataset | floor hard-recall | mutual-kNN (k=16) hard-recall |
|---|---|---|
| seoul_bull (17 imgs, perspective) | 86% | 100% |
| seattle_backyard (26) | 54% | 99% |
| kerry_park (48, back-to-back fisheye) | 52% | 95% |
| dino_dog (85, object-centric, dense texture) | 30% | 68% |

The floor leaves roughly **half** the hard-pair matches on the table on the
harder collections; mutual-kNN roughly **doubles** hard-pair recall everywhere.
`dino_dog` is the ceiling — a crowded descriptor space pushes true matches past
rank-16, so it benefits from a larger `k` — but the *relative* improvement is the
same on every dataset.

### The recall gain becomes reconstruction density

End-to-end (mutual-kNN candidates -> geometric verification -> global solve),
median reconstructed points over repeated solves:

| dataset | exhaustive | cluster (floor) | mutual-kNN (k=12) |
|---|---|---|---|
| seoul_bull | 732 | 853 | **948** (+11% vs floor) |
| kerry_park | 402 | 349 | **420** (+20%) |
| seattle_backyard | 3205 | 3327 | **4800** (+44%) |

Mutual-kNN produces the densest reconstruction on every dataset, matching or
beating exhaustive (it beats exhaustive because COLMAP's exhaustive matcher
applies its own ratio test; mutual-kNN keeps more true matches) at clustering
cost.

### It is a density / completeness gain, not an accuracy gain

`sfm compare` between the exhaustive and mutual-kNN reconstructions shows the
**cameras are essentially unchanged** (on well-conditioned scenes the poses
agree to < 0.1% of scene scale); the mutual-kNN cloud is the exhaustive cloud
**plus** the recovered wide-baseline points. Those extra points are **noisier** —
on seattle the median reprojection error rises 0.30 -> 0.41 px and the 90th
percentile 0.70 -> 0.99 px. So mutual-kNN trades a measurable increase in point
noise for ~40% more points and completeness. Treat it as a **completeness**
knob.

### What does *not* work

- **Lowe's ratio test and mutual-top-1** collapse hard-pair recall (to 2–8%).
  The classic two-view ratio test assumes a feature has *one* correct match;
  a multi-view track has several near-equal correct neighbours (the other
  views), so the ratio test wrongly rejects them. Ratio/top-1 are precision
  filters tuned for the opposite regime and must not be applied here.
- **Connected components / transitive closure** of the mutual graph over-merges:
  a few ambiguous bridges chain large swaths of the corpus into one giant
  component (76% of the corpus on kerry_park at k=6), exploding the candidate
  count. Mutual-kNN therefore emits the **edges directly** — it never takes the
  transitive closure.

### The triangle filter (optional precision)

A real 3-D point seen in >= 3 images has observations that are pairwise mutually
near, so they form a **triangle (3-clique)** in the mutual graph. A spurious
match needs a third descriptor that is *independently* mutually near to *both*
endpoints — unlikely for a random collision. The triangle filter keeps an edge
only if it closes at least `triangle_min` triangles. It is a geometry-free,
multi-view-consistency pre-filter: effectively "keep this match only if it is
part of a >= 3-view mini-track."

On kerry_park (k=12), it is a clean recall/precision dial:

| `triangle_min` | candidate matches | hard-recall | spurious |
|---|---|---|---|
| 0 (off) | 229k | 87% | 81% |
| 1 | 150k | 82% | 72% |
| 2 | 97k | 75% | 61% |
| 3 | 68k | 66% | 50% |
| *(floor, for comparison)* | *53k* | *52%* | *36%* |

At `triangle_min=3` it reaches the floor's verification budget (candidate count
and spurious rate) but with strictly better hard-recall (66% vs 52%); at
`triangle_min=1` it keeps most of the recall while shedding ~35% of the
candidates. The third triangle corner is always a genuine third image (the k-NN
excludes same-image hits), so a surviving wide edge is one corroborated by an
intermediate view — exactly the ones most likely to be real.

> Note: do **not** use triangles to *add* edges (closing 2-hop paths through a
> shared mutual neighbour). That recovers a little more recall but reintroduces
> the connected-components over-merge — on kerry_park it pushed the candidate
> count to ~1M at 96% spurious. The triangle step is a *filter* only.

## Algorithm

Given the `(N, 128)` uint8 descriptor corpus (every image's SIFT descriptors
concatenated image by image) and a CSR `image_starts` array of length
`n_images + 1` mapping rows to images:

1. **Index + shared k-NN query.** Build one `KdForest` over the corpus and run a
   single batched k-NN query (`search_batch_with_distances`). Because the forest
   returns nearest neighbours over the whole corpus — including self and
   same-image hits — query a generous width `2k + 1` so that, after dropping
   self and same-image rows, each descriptor still has its `k` nearest genuine
   cross-image neighbours.
2. **Mutual edges.** For each descriptor `a` and each of its `k` cross-image
   neighbours `b` with `a < b`, keep `a-b` iff `a` is also among `b`'s `k`
   neighbours. Membership is a short linear scan of `b`'s length-`k` row.
3. **Triangle filter (optional).** If `triangle_min > 0`, build the undirected
   adjacency of the surviving edges and keep `a-b` only if `|adj(a) ∩ adj(b)| >=
   triangle_min`.
4. **Expand to per-image-pair matches.** Because the corpus is concatenated
   image by image, for an edge `a < b` in different images `image_of[a] <
   image_of[b]`, so the image pair is already ordered. Emit, bucketed by image
   pair and sorted, the same four parallel arrays the `.matches` writer wants:
   `image_index_pairs (P,2)`, `match_counts (P,)`, `match_feature_indexes (M,2)`,
   `match_descriptor_distances (M,)`. This is the identical shape
   `clusters_to_pair_matches` returns, so the rest of the matching pipeline
   (DB population, geometric verification, `.matches` writing) is shared.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `k` | 12 | Nearest cross-image neighbours kept per descriptor. Larger `k` recovers more wide-baseline matches (higher recall) at more candidates to verify. 12–16 is the sweet spot on the bundled datasets; very dense / object-centric collections benefit from higher `k`. |
| `triangle_min` | 0 | Triangle filter: keep an edge only if it closes >= this many triangles. 0 disables. 1–2 sheds spurious candidates with little recall loss; 3 reaches floor-like verification cost at better recall. |
| `preset` (forest) | `accurate` | Kd-forest build + per-query budget. `accurate` keeps the k-NN faithful enough that mutuality is meaningful; `fast`/`balanced` trade recall for speed. |

## Cost Analysis

One `KdForest` build and one batched k-NN query, both `O(N log N)` and shared
across the whole collection — the same cost class as the background-floor
matcher (which already builds the forest and queries it). The mutual-edge and
triangle passes are `O(N · k)` and `O(E · k)`. This is **not** the `O(n_images^2)`
all-pairs descriptor matching of exhaustive: the per-image-pair structure falls
out of the global k-NN. The extra cost relative to the floor matcher is paid
**downstream**, in geometric verification, because mutual-kNN hands it more
candidate matches; the triangle filter exists to bound that.

## Relationship to the Existing Pipeline

Mutual-kNN is a sibling of the background-floor matcher: both are corpus-level
matchers that emit per-image-pair matches into the shared `.matches` pipeline
(`PairMatches` / `PairArrays` -> COLMAP DB -> `verify_matches` -> `.matches`
writer). They differ only in how candidates are selected from the shared k-NN
query — adaptive radius + cluster partition (floor) vs rank-bounded mutual edges
(this matcher). The floor optimizes precision (fewer, cleaner candidates);
mutual-kNN optimizes recall/completeness and delegates precision to RANSAC.

## Limitations

- **Verification cost.** The high candidate count (3–6x the floor's) means more
  two-view geometry estimations. Mitigate with the triangle filter or a smaller
  `k`.
- **Point noise.** The recovered wide-baseline points are noisier (see above).
  For accuracy-critical use prefer the floor matcher or a stricter
  `triangle_min`.
- **`k` is collection-dependent.** Crowded descriptor spaces (large, repetitive,
  object-centric collections) push true matches past a fixed `k`; recall there is
  `k`-limited.
- **Approximate k-NN.** The forest is approximate, so the realized neighbour sets
  (and thus recall) depend on the forest preset; `accurate` is the default for
  that reason.

## Implementation

### Layer 1 — Rust core (`sfmtool-core`)

`crates/sfmtool-core/src/mutual_knn/mod.rs`.

```rust
pub struct MutualKnnParams { pub k: usize, pub triangle_min: usize, pub forest: KdForestParams }

pub fn mutual_knn_matches(
    descriptors: ArrayView2<'_, u8>,
    image_starts: &[u32],
    params: &MutualKnnParams,
) -> Result<PairMatches, ClusterMatchError>;
```

Reuses `PairMatches` and `ClusterMatchError` from `crate::cluster_match` (the
output shape and input-validation errors are identical). Parallelised with
rayon (the k-NN query, the mutual-edge pass, the triangle filter, and the
expansion). Deterministic given a fixed forest seed: edges are produced in a
fixed order and the final tuples are sorted by `(img_i, img_j, feat_i, feat_j)`.
Unit tests in `crates/sfmtool-core/src/mutual_knn/tests.rs`.

### Layer 2 — PyO3 binding (`sfmtool-py`)

`crates/sfmtool-py/src/py_mutual_knn.rs` exposes `mutual_knn_matches(descriptors,
image_starts, k=12, triangle_min=0, preset=None, ...)` returning the four arrays
`(image_index_pairs, match_counts, match_feature_indexes,
match_descriptor_distances)` — the same tuple `clusters_to_pair_matches`
returns.

### Layer 3 — Python (`sfmtool`)

- `src/sfmtool/feature_match/_mutual_knn_matching.py`: `mutual_knn_match(...)`
  builds the corpus from the `.sift` files and calls the binding, returning a
  `PairArrays` (reused from `_cluster_matching`).
- `src/sfmtool/feature_match/_run.py`: `_run_mutual_knn_matching(...)` runs the
  matcher and writes the result to the DB via the shared
  `_write_pairs_and_verify(...)` helper (factored out of the cluster path); the
  `"mutual-knn"` method is wired into `_run_matching` and recorded in the
  `.matches` metadata (`mode="mutual-knn"`, `k`, `triangle_min`, `preset`).

### Layer 4 — CLI

`sfm match --mutual-knn [--mutual-knn-k K] [--mutual-knn-triangle T]
[--mutual-knn-preset P]` (`src/sfmtool/_commands/match.py`). See
`specs/cli/match-command.md`.
