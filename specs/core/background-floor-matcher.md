# Background-Floor Matcher — Production Implementation Spec

This spec describes the production form of the **background-floor** track-cluster
matcher (the `bgfloor` mode of the POC `experiments/exp05_cluster_match.py`),
across three layers: a Rust matcher in `sfmtool-core`, its PyO3 binding, and the
`sfm match` CLI wiring. It is written so a fresh implementation can be done
against it without reading the POC, though
[the POC reference table](#appendix-poc-reference) maps each piece back to the
experiment code.

The **algorithm and its empirical justification** are specified separately in
[`track-cluster-matching.md`](track-cluster-matching.md), section *Membership by
per-point background floor (recommended)*. This document does not re-derive the
method; it pins down the API and data flow. Read that section first for the
"why".

## What the matcher does (one paragraph)

Given the SIFT descriptors of every image in a set, concatenate them into one
corpus, build a randomized kd-tree forest over it, and query each descriptor's
`k` nearest neighbours. For each descriptor, set a *per-point* radius from its own
background floor — `alpha ×` the median neighbour distance from rank `b0` onward —
and keep its cross-image neighbours within that radius. Reconcile the resulting
directed edges to one match per feature per image pair, and emit a `.matches`
file (no two-view geometry) for the existing `sfm solve` / `sfm to-colmap-db`
consumers. There is no global threshold and no pose/intrinsics input.

## Distance space (read this first)

All distances in this matcher are **Euclidean L2** (square-rooted), not squared.
This matters because the tuned defaults `alpha = 0.8`, `b0 = 8` were fit in L2
space (via Python `KdForest.query`, which returns L2).

The core `KdForest::search_batch_with_distances` returns **squared** L2
(`dist_sq`). **The matcher must take the square root of those distances before
computing the background-floor median and the radius test.** Distances written to
`match_descriptor_distances` are likewise L2, matching every other matcher's
`.matches` output.

---

## Layer 1 — Rust core (`sfmtool-core`)

### Location

New module `crates/sfmtool-core/src/cluster_match/` with `mod.rs`; declare
`pub mod cluster_match;` in `crates/sfmtool-core/src/lib.rs`. Optionally re-export
the public types from `lib.rs` alongside the other `pub use` lines.

### Public types

```rust
use ndarray::{Array1, Array2, ArrayView2};
use crate::kdforest::KdForestParams;

/// Tuning for the background-floor matcher. `Default` is the production config.
#[derive(Clone, Debug)]
pub struct BackgroundFloorParams {
    /// Neighbours fetched per descriptor, including self (column 0). Default 49.
    pub k: usize,
    /// Radius multiplier: keep neighbours within `alpha * B_i`. Default 0.8.
    pub alpha: f32,
    /// First neighbour rank counted as background when estimating `B_i`
    /// (`B_i = median(dist[i, b0..k])`). Default 8.
    pub b0: usize,
    /// Index build + per-query search budget. Default `KdForestParams::accurate()`.
    pub forest: KdForestParams,
}

impl Default for BackgroundFloorParams {
    fn default() -> Self {
        Self { k: 49, alpha: 0.8, b0: 8, forest: KdForestParams::accurate() }
    }
}

/// Cross-image matches, in the parallel-array form the `.matches` writer wants.
pub struct PairMatches {
    /// `(P, 2)` image-index pairs, each `[i, j]` with `i < j`, sorted ascending
    /// by `(i, j)`.
    pub image_index_pairs: Array2<u32>,
    /// `(P,)` number of matches in each pair; `sum == M`. Aligned to
    /// `image_index_pairs`.
    pub match_counts: Array1<u32>,
    /// `(M, 2)` feature-index pairs `[feat_i, feat_j]`, grouped by pair in the
    /// same order as `image_index_pairs`. `feat_i` indexes image `i`'s `.sift`
    /// rows, `feat_j` indexes image `j`'s.
    pub match_feature_indexes: Array2<u32>,
    /// `(M,)` Euclidean L2 descriptor distance per match, aligned to
    /// `match_feature_indexes`.
    pub match_descriptor_distances: Array1<f32>,
}

#[derive(Debug, thiserror::Error)]  // or a hand-rolled enum, matching crate style
pub enum ClusterMatchError {
    #[error("descriptor corpus is empty")]
    EmptyCorpus,
    #[error("k ({k}) must be > b0 ({b0}) so the background median has samples")]
    BadK { k: usize, b0: usize },
    #[error("image_starts must be non-decreasing, start at 0, and end at N ({n})")]
    BadOffsets { n: usize },
}
```

> Use whichever error idiom the crate already uses — check neighbouring modules
> (e.g. `reconstruction`, `camera_intrinsics`) and match it; `thiserror` above is
> illustrative.

### Public function

```rust
/// Background-floor track-cluster matcher (spec: track-cluster-matching.md,
/// "Membership by per-point background floor").
///
/// `descriptors` is the `(N, D)` corpus of every image's SIFT descriptors,
/// concatenated image by image (uint8, D = 128). `image_starts` is a CSR-style
/// offset array of length `n_images + 1`: image `i` owns rows
/// `image_starts[i] .. image_starts[i+1]`, and row `r` of that image has
/// feature index `r - image_starts[i]`. Returns one-to-one-per-image-pair
/// cross-image matches.
pub fn background_floor_match(
    descriptors: ArrayView2<'_, u8>,
    image_starts: &[u32],
    params: &BackgroundFloorParams,
) -> Result<PairMatches, ClusterMatchError>;
```

### Algorithm (exact)

Let `N = descriptors.nrows()`, `D = descriptors.ncols()` (128),
`n_images = image_starts.len() - 1`.

1. **Validate.** `N > 0`; `params.k > params.b0`; `image_starts[0] == 0`,
   non-decreasing, `image_starts[n_images] == N`. Else return the matching
   `ClusterMatchError`.

2. **Row → (image, feature) maps.** From `image_starts`, build `image_of[r]`
   (`u32`, the owning image) and `feature_of[r]` (`u32`, `r - image_starts[image]`)
   for every row `r`. (A binary search over `image_starts`, or a single linear
   fill, both fine.)

3. **Build the forest.** `KdForest::build` over the flat `descriptors` slice
   (row-major `N*D` `u8`), `dim = D`, with `params.forest`. The corpus must be
   contiguous; if `descriptors` is not standard layout, copy to a `Vec<u8>` first.

4. **Query.** `let (idx, dist_sq) = forest.search_batch_with_distances(corpus, N,
   params.k, params.forest.max_leaf_checks, None);` — flat `N*k` arrays, each row
   sorted ascending with column 0 = self at distance 0. Unfound slots are
   `u32::MAX` / `f32::INFINITY` (only possible if an image has fewer than `k`
   features total across the corpus — corpus is global, so rare).

5. **L2 distances.** Materialise `dist[r*k + c] = dist_sq[...].sqrt()`.

6. **Per-point radius.** For each row `i`: `B_i = median(dist[i, b0..k])` (the
   `k - b0` farthest of the fetched neighbours; with the defaults, 41 samples,
   none of them self). `radius_i = alpha * B_i`. Use a simple median (sort the
   slice or `select_nth`; even count → average the two middle values, matching
   NumPy's `np.median`).

7. **Candidate edges.** For each row `i` and each neighbour column `c` in
   `0..k` with `j = idx[i, c]`: keep the directed edge `i → j` iff
   `j != u32::MAX`, `j != i`, `image_of[i] != image_of[j]`, and
   `dist[i, c] <= radius_i`. Record `(i, j, dist[i,c])`. (Self at column 0 is
   dropped by `j != i`.)

8. **Canonicalise.** For each kept edge map to
   `(img_lo, img_hi, feat_lo, feat_hi, d)` where `img_lo < img_hi`: if
   `image_of[i] < image_of[j]` then `(image_of[i], image_of[j], feature_of[i],
   feature_of[j])` else the swap. `d` is the L2 distance (symmetric, so order
   does not change it).

9. **One-to-one per image pair (two-pass min).** Within each image pair, keep at
   most one match per `feat_lo` and per `feat_hi`, preferring smaller `d`:
   - **Pass 1:** group by `(pair, feat_lo)`; keep only the minimum-`d` edge in
     each group.
   - **Pass 2:** on the survivors, group by `(pair, feat_hi)`; keep only the
     minimum-`d` edge in each group.

   This is the POC's two-pass `keep_min_by_key` and guarantees one-to-one per
   image pair (each low-feature and each high-feature appears at most once per
   pair). It is deliberately cheaper than global mutual-best; geometric
   verification downstream enforces strict geometric one-to-one. Implement with a
   stable sort by the composite key then a "first per key after sorting by `d`"
   reduction, or an equivalent hash reduction.

10. **Assemble.** Sort the surviving edges by `(img_lo, img_hi)`. Produce:
    `image_index_pairs` = the distinct sorted pairs; `match_counts` = per-pair
    edge counts; `match_feature_indexes` = `[feat_lo, feat_hi]` rows grouped by
    pair in that order; `match_descriptor_distances` = the aligned `d` values.

### Parallelism

Steps 6–8 are embarrassingly parallel over rows — use `rayon` (`par_iter` /
`par_chunks`) as elsewhere in the crate. Step 4's `search_batch_with_distances`
is already internally parallel. The step-9 reduction can use a parallel sort
(`rayon`'s `par_sort_unstable_by`). Keep memory bounded: the candidate-edge set
is `≤ N*(k-1)`; for dino (`N ≈ 600k`, `k = 49`) that is ~29M edges before dedup —
fine as flat `Vec`s of primitives, but do not build per-edge structs with heap
fields.

### Determinism

Given a fixed `params.forest.seed`, the matcher is deterministic. Break `d` ties
in steps 6/9 deterministically (e.g. by neighbour index) so output ordering is
stable across runs and platforms.

### Tests (`crates/sfmtool-core/src/cluster_match/tests.rs`)

- **Synthetic clusters.** Build a tiny corpus: a few "points" each with one
  descriptor in 3–4 images (tight intra-cluster distance) plus scattered
  background. Assert the matcher recovers the cross-image pairs of each planted
  point and emits no within-image pairs.
- **One-to-one.** Assert that within every returned image pair, no `feat_lo` and
  no `feat_hi` repeats.
- **Distances are L2.** Assert a returned distance equals the true L2 between the
  two descriptors (within forest approximation / float tolerance), not squared.
- **Validation.** `k <= b0`, empty corpus, and malformed `image_starts` return
  the right `ClusterMatchError`.
- **Determinism.** Two runs with the same seed produce byte-identical arrays.

Run `pixi run cargo test -p sfmtool-core cluster_match` and
`pixi run cargo clippy --workspace` / `pixi run cargo fmt`.

---

## Layer 2 — PyO3 binding (`sfmtool-py`)

### Location

New file `crates/sfmtool-py/src/py_cluster_match.rs`; `mod py_cluster_match;` and
registration in `crates/sfmtool-py/src/lib.rs`’s `#[pymodule]` via
`m.add_function(wrap_pyfunction!(py_cluster_match::background_floor_match, m)?)?;`.

### Function

Mirror the `KdForest` binding conventions (`py_kdforest.rs`): validate the uint8
dtype, make the corpus contiguous, release the GIL around the heavy call, and
return numpy arrays via `PyArray1::from_vec(...).reshape(...).into_any().unbind()`.

```rust
/// Background-floor track-cluster matcher.
///
/// Args:
///     descriptors: (N, 128) uint8 corpus, every image's SIFT descriptors
///         concatenated image by image.
///     image_starts: (n_images + 1,) uint32 CSR offsets; image i owns rows
///         image_starts[i]:image_starts[i+1].
///     k, alpha, b0: background-floor parameters (defaults 49, 0.8, 8).
///     preset / num_trees / leaf_size / max_leaf_checks / seed: forest config,
///         same meaning as KdForest.
///
/// Returns:
///     (image_index_pairs (P,2) uint32, match_counts (P,) uint32,
///      match_feature_indexes (M,2) uint32, match_descriptor_distances (M,) float32)
#[pyfunction]
#[pyo3(signature = (descriptors, image_starts, k=49, alpha=0.8, b0=8,
                    preset=None, num_trees=None, leaf_size=None,
                    max_leaf_checks=None, seed=None))]
pub fn background_floor_match<'py>(
    py: Python<'py>,
    descriptors: &Bound<'py, PyAny>,
    image_starts: PyReadonlyArray1<'py, u32>,
    k: usize,
    alpha: f32,
    b0: usize,
    preset: Option<&str>,
    num_trees: Option<usize>,
    leaf_size: Option<usize>,
    max_leaf_checks: Option<usize>,
    seed: Option<u64>,
) -> PyResult<(Py<PyAny>, Py<PyAny>, Py<PyAny>, Py<PyAny>)>;
```

Build `KdForestParams` from `preset` + overrides exactly as `PyKdForest::new`
does (reuse that resolution logic — consider lifting it into a shared helper).
Map `ClusterMatchError` to `PyValueError::new_err(...)`. Validate
`descriptors.ndim == 2`, `ncols == 128`, dtype `uint8`, and
`image_starts.len() == n_images + 1` with a clear message.

### Python package surface

`KdForest` is re-exported as `sfmtool.KdForest` (see `src/sfmtool/__init__.py`).
Re-export the new function the same way so callers can `from sfmtool import
background_floor_match` (it lives in `sfmtool._sfmtool`).

### Rebuild + tests

After the Rust edits, **`pixi run maturin develop --release`** (the `.so` does
not auto-rebuild). Add `tests/rust_bindings/test_cluster_match_rust_bindings.py`:

- A tiny hand-built corpus (numpy) with a couple of planted cross-image points;
  assert the returned arrays have the documented shapes/dtypes, that pairs are
  sorted and `i < j`, that `match_counts.sum() == len(match_feature_indexes)`,
  and that planted matches appear.
- Dtype/shape errors raise `ValueError`/`TypeError`.

---

## Layer 3 — Python matcher layer + CLI

### Matcher module

New `src/sfmtool/feature_match/_cluster_matching.py`, mirroring
`_flow_matching.py`:

```python
def cluster_match(
    image_paths: list[Path],
    sift_paths: list[Path],
    *,
    k: int = 49,
    alpha: float = 0.8,
    b0: int = 8,
    preset: str = "accurate",
    max_feature_count: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Run the background-floor matcher over every image's SIFT descriptors.

    Loads each image's descriptors (capped at max_feature_count to match the
    feature indices used downstream), concatenates them into one (N, 128) uint8
    corpus with a CSR image_starts array, and calls
    sfmtool.background_floor_match. Returns the four parallel arrays
    (image_index_pairs, match_counts, match_feature_indexes,
    match_descriptor_distances).
    """
```

Load descriptors with `SiftReader(sift_path).read_descriptors(count=...)` (see
`src/sfmtool/sift/file.py`). Build `image_starts` as the cumulative feature
counts. Feature indices in the result are `.sift` row indices (capped to the
first `max_feature_count` when set), consistent with how `to-colmap-db` loads
keypoints.

### Writing the `.matches` file

Reuse the existing assembly used by `feature_match/_run.py` (`_run_matching`):
build the `matches_data` dict with `has_two_view_geometries = False`, per-image
`feature_tool_hashes` / `sift_content_hashes` / `feature_counts` from
`read_sift_metadata`, workspace metadata pulled from `.sfm-workspace.json`, and
the four pair arrays; then `from sfmtool._sfmtool import write_matches;
write_matches(out, matches_data)`. Set `matching_method = "cluster"`,
`matching_tool = "sfmtool"`, `matching_tool_version` = the package version, and
record the parameters in `matching_options` (`{"mode": "background-floor", "k":
k, "alpha": alpha, "b0": b0, "preset": preset}`).

> Factor the dict assembly so the cluster path and the existing paths share it
> rather than duplicating the metadata plumbing. The POC's
> `assemble_matches_dict` is the shape to reproduce, but use the package's real
> tool name/version, not the POC's `"sfmtool-cluster-poc"`.

Add an orchestrator in `_run.py`, e.g. `_run_cluster_matching(image_paths,
output_path, *, k, alpha, b0, preset, max_feature_count, workspace_dir)`, that
resolves `.sift` paths via `image_files_to_sift_files`, calls `cluster_match`,
assembles the dict, writes it, and returns the output path. Default output to
`workspace/matches/<timestamp>-cluster.matches` (cluster matches never carry
two-view geometry, so they go under `matches/`, not `tvg-matches/`).

### CLI: extend `sfm match`

Add a fourth matching **method** to `src/sfmtool/_commands/match.py`, mutually
exclusive with `--exhaustive` / `--sequential` / `--flow`:

```
--cluster                 use the background-floor track-cluster matcher
--cluster-alpha FLOAT     background-floor radius multiplier (default 0.8)
--cluster-b0 INT          first neighbour rank counted as background (default 8)
--cluster-k INT           neighbours fetched per descriptor incl. self (default 49)
--cluster-preset CHOICE   forest preset: accurate|balanced|fast (default accurate)
```

`--cluster` dispatches to `_run_cluster_matching` over the resolved image set and
honours `--max-features`, `--output`, and `--range` (range restricts the image
set the corpus is built from). Update the method-selection / mutual-exclusion
validation and the help text. Keep it in the **Image Feature** category (already
where `match` is registered in `cli.py`).

**Camera model / camera_config.** The cluster matcher uses no intrinsics or
poses, so `--camera-model` is not applicable to `--cluster`; reject the
combination with a `click.UsageError` (it only means something for the
registered-image descriptor matcher). The existing
`_check_camera_model_conflict` (`src/sfmtool/_camera_setup.py`) still applies to
the other methods unchanged.

### Spec

Add a short `## Cluster matching` section to `specs/cli/match-command.md`
describing `--cluster` and its options, and link back to this file and to
`track-cluster-matching.md`.

### Downstream (no new code)

The emitted `.matches` flows through the existing consumers unchanged:

```bash
pixi run sfm match --cluster images -o matches/cluster.matches
pixi run sfm to-colmap-db matches/cluster.matches --out-db colmap.db   # then verify + map
pixi run sfm solve -i matches/cluster.matches                          # incremental SfM
```

`sfm to-colmap-db` / `sfm solve` already turn `has_two_view_geometries = false`
matches into a COLMAP DB and run COLMAP's own geometric verification per pair
(`_setup_for_sfm_from_matches` / `_write_matches_to_db` in `_colmap_db.py`). A
`sfm solve --cluster` shortcut (match then map in one call) is a reasonable
follow-on but is **out of scope** here; the `match → .matches → solve` path is
the deliverable.

### Tests

- **Unit**: `tests/test_cluster_matching.py` — small synthetic descriptor sets
  through `cluster_match`, asserting shapes and one-to-one-per-pair.
- **Integration**: using the `isolated_seoul_bull_17_images` fixture (see
  `tests/conftest.py`), run `sfm match --cluster` and assert a `.matches` file is
  produced with the expected pair/match counts > 0 and `has_two_view_geometries`
  False; optionally feed it to `to-colmap-db` and assert a DB is created.
- `pixi run fmt && pixi run check` for the Python changes.

---

## Defaults (single source of truth)

| Parameter | Default      | Layer(s)            | Meaning                                            |
| --------- | ------------ | ------------------- | -------------------------------------------------- |
| `k`       | 49           | core/py/cli         | neighbours fetched incl. self                      |
| `alpha`   | 0.8          | core/py/cli         | keep neighbours within `alpha · B_i`               |
| `b0`      | 8            | core/py/cli         | first neighbour rank counted as background         |
| `preset`  | `accurate`   | core (`KdForestParams`)/py/cli | forest build + search budget            |
| distance  | Euclidean L2 | all                 | sqrt of the forest's squared distances             |
| TVG       | none         | py/cli              | `has_two_view_geometries = false`                  |

These are the tuned values from `track-cluster-matching.md`; do not change them
without re-running the exp20–23 bench and the end-to-end reconstructions.

## Implementation order (suggested)

1. Rust `cluster_match` module + unit tests → `cargo test`/`clippy`/`fmt`.
2. PyO3 binding + `maturin develop --release` + `tests/rust_bindings/...`.
3. `_cluster_matching.py` + `_run.py` orchestrator + Python unit test.
4. `sfm match --cluster` wiring + integration test + `specs/cli` note.
5. `fmt && check`, full `pixi run test`, update this file's status if anything
   diverged.

## Appendix: POC reference

Every piece exists in `experiments/exp05_cluster_match.py`; the production job is
to move it into the crate/bindings/CLI with real metadata and tests.

| Production element                         | POC source                                             |
| ------------------------------------------ | ------------------------------------------------------ |
| corpus + `image_starts` + (image,feature)  | `sfm_descriptors.load_descriptor_bank` (concat blocks) |
| `k = 49` query                             | `BG_K = 49`; `forest.query(..., k=query_k)`            |
| per-point radius `alpha * median(dist[b0:])`| `build_bgfloor_matches_arrays` (`radius = alpha*np.median(dst[:, b0:], axis=1)`) |
| edge keep test                             | `keep = (dd <= radius[i_rep]) & (i_rep!=j) & (img[i_rep]!=img[j])` |
| two-pass one-to-one dedup                  | `_edges_to_pair_arrays` + `keep_min_by_key`            |
| pair arrays                                | `_edges_to_pair_arrays` return                         |
| matches dict + write                       | `assemble_matches_dict` + `write_matches`              |
| end-to-end driver (reference behaviour)    | `experiments/run_bgfloor.sh`                           |

Note the POC writes the forest distance directly as the match distance; since the
binding returns L2, that is already L2 — keep it L2 in production too.
