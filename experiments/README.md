# Feature-matching experiments (POC)

Standalone proof-of-concept code for exploring feature-matching ideas against
the checked-in test datasets. **Not** wired into the `sfmtool` package — these
scripts read existing reconstructions and `.sift` files directly.

Everything runs in the dedicated `experiments` pixi environment (adds
scikit-learn + matplotlib):

```bash
pixi run -e experiments python experiments/exp01_descriptor_kdtree_clusters.py \
    seoul_bull_ws/sfmr/*.sfmr
```

## Setup used

All four datasets are bootstrapped from `scripts/init_dataset_*.sh`, then solved:

| Dataset            | Workspace             | Imgs | Solve       | Points | Character                          |
| ------------------ | --------------------- | ---- | ----------- | ------ | ---------------------------------- |
| seoul_bull_sculpt. | `seoul_bull_ws/`      | 17   | incremental | 1,080  | close-range textured, small bases  |
| seattle_backyard   | `seattle_backyard_ws/`| 26   | global      | 521    | outdoor, moderate baselines        |
| kerry_park         | `kerry_park_ws/`      | 48   | global      | 1,128  | 2-camera fisheye rig, distortion   |
| dino_dog_toy       | `dino_dog_toy_ws/`    | 85   | global      | 5,312  | turntable, wide baselines          |

(Workspaces are gitignored under `/*_ws/`; re-create with the init + solve
scripts.)

## `sfm_descriptors.py`

Shared loader. `load_descriptor_bank(sfmr_path)` returns a `DescriptorBank`:
every SIFT descriptor from every image concatenated into one `(N, 128)` array,
each row labelled with the solve's 3D-point id of the track it belongs to
(`-1` if the solve never used it). The solve's tracks are the *ground truth*.

## exp01 — do descriptor-space clusters look like the solve's tracks?

The original idea: throw all descriptors into one ANN/KD-tree and see whether
the structure of descriptor space recovers the feature tracks the solve built.
Four probes:

- **A. Separation** — within-track vs random cross-image distances.
- **B. k-NN recovery** — are a feature's co-track members among its NN?
- **C. DBSCAN vs tracks** — density clustering scored against track labels.
- **D. Mutual-NN components** — cross-image mutual NN → connected components.

### Findings (2026-06-06, all four datasets)

**The structure is there, and its strength tracks scene difficulty.**

| Dataset          | Descriptors | In-track | Top-1 NN = match | recall@5 | Mutual-NN purity | Best DBSCAN (homog/compl) |
| ---------------- | ----------- | -------- | ---------------- | -------- | ---------------- | ------------------------- |
| seoul_bull       | 37.8k       | 10.6%    | **98.7%**        | 0.89     | 0.987            | 0.86 / 0.96 (eps≈107)     |
| seattle_backyard | 79.9k       | 3.5%     | **93.4%**        | 0.71     | 0.945            | 0.90 / 0.96 (eps≈135)     |
| kerry_park       | 81.0k       | 5.8%     | **81.7%**        | 0.74     | 0.856            | 0.91 / 0.96 (eps≈109)     |
| dino_dog_toy     | 336.5k      | 11.0%    | **81.5%**        | 0.46     | 0.858            | 0.90 / 0.98 (eps≈104)     |

- **A:** within-track distances sit far below random cross-image pairs in every
  dataset (medians ~105–135 vs ~270–380); 88–99% of true co-track distances
  fall below the random 5th percentile. Clear signal everywhere.
- **B (the headline):** the single nearest cross-image descriptor is the solve's
  matched feature **82–99%** of the time — an excellent, nearly free *pairwise*
  match signal straight off the index.
- **C:** DBSCAN works only in a narrow `eps` band around the within-track median
  (homogeneity ~0.86–0.91, completeness ~0.96), and collapses into a few
  mega-clusters just above it (and OOMs at large eps on the big set). `eps` is
  make-or-break.
- **D:** mutual-NN pairs are **86–99% pure** with **100%** one-feature-per-image
  — exactly the structure real tracks have. *But* pure mutual-1NN only ever
  yields **pairs**, so the solve's 3–9-image tracks come out fragmented into
  disjoint pairs rather than aggregated into full tracks.

**What the dataset differences tell us:**

- **Scene difficulty degrades the signal predictably.** seoul_bull (close-range,
  textured, small baselines) is near-perfect; **kerry_park's fisheye distortion**
  and **dino's wide-baseline turntable views** both drop top-1 to ~81%.
- **Wide baselines wreck k-NN *recall* specifically** (dino recall@5 = 0.46): a
  track's many wide-baseline observations scatter across descriptor space, so
  they are *not* all mutual near-neighbours — which is exactly why mutual-1NN
  fragments tracks, and worse the wider the baselines.
- Only **3.5–11%** of all descriptors end up in any solve track — descriptor
  space is dominated by features the solve discarded, so clustering over *all*
  descriptors over-segments (far more clusters than tracks).

**Implementation notes:**

- A literal KD-tree is pathological at 128-D (sklearn's `kd_tree` hung for 15+
  CPU-min on the 80k set). These analysis probes compute neighbours exactly via a
  chunked BLAS GEMM (`knn.py`) as the oracle; the 336k×336k exact mutual-NN takes
  ~27 min. The matcher (`exp05`) uses the in-tree approximate index
  (`sfmtool.KdForest`) instead, which is the production path.
- DBSCAN doesn't scale past ~100k points; for large sets the probe subsamples
  the corpus (keeping every in-track descriptor so scoring is unaffected) and
  caps the eps sweep below the always-collapses regime.

**Where this points next:** descriptor NN gives near-perfect *pairs* cheaply;
the open problem is **aggregating pairs into full multi-image tracks**
(transitive closure / graph clustering with the one-per-image constraint),
then feeding the resulting matches into `sfm solve` to measure the downstream
reconstruction vs the baselines above.

## exp02 — neighbour-distance structure (is there a threshold?)

For a dino sample, computes exact 6-NN with true distances and asks whether the
nearest distance `d1` or the `d1/d5` ratio separates real-track descriptors from
background. Plots are written to `out/`.

**Findings.** Both are *rejectors, not selectors*:

- `d1` is **bimodal**: a "has-a-close-neighbour" mode (~45–55) that overlaps
  in-track, and an "isolated" mode (~100–110) that is background-only.
- the `d1/d5` ratio has a sharp **high mode (~0.9–0.97)** that is background-only
  ("my 5 nearest are all equally far" = isolated).
- So `ratio > ~0.85` or `d1` above the antimode cheaply drops ~40% of
  descriptors as isolated background. **But** a large population of background
  descriptors have a genuinely close neighbour (low `d1`, mid ratio) and are
  *indistinguishable from real track members by distance alone* — these are the
  repeated-structure / failed-verification matches that only geometry removes.

## exp03 — bounded-radius clustering vs the solve's tracks

Pipeline: for every descriptor take its ≤16 nearest, optionally drop isolated
points (the exp02 prefilter), keep cross-image edges within radius `T`, take
connected components as candidate tracks, and score them against the solve's
tracks over a `T` sweep. The 17-NN is computed once (exact) and cached so the
sweep is cheap.

**Findings (dino).** A real candidate-track generator:

- Sweet spot **T≈80 + prefilter**: purity **0.98**, track completeness **0.87**,
  **1.6** fragments/track, 95% image-unique. A clear step beyond exp01's
  mutual-NN *pairs* — the radius graph aggregates into multi-image components.
- The **isolated-point prefilter is a strict Pareto win** (better recovery *and*
  purity at the same `T`).
- There is a **sharp purity/recovery knee**: past T≈100 purity collapses (false
  merges through the overlap-region background). Crossing it needs geometry.
- Residual gaps are all geometry-shaped: ~13% of members never join a component
  (wide-baseline scatter), ~5% of components mix two features from one image,
  and safe merging past the knee needs verification.

## exp04 — deriving the radius from the data (no labels)

exp03's `T` was hand-picked. But `d1` is bimodal, so the **antimode is a
label-free threshold**. We pick it two ways — Otsu and a 2-component Gaussian
mixture — then validate the resulting clustering against the solve's tracks.

**Findings — the threshold must be, and can be, derived from data.** The
absolute distance scale varies ~2.5× across datasets, so a fixed radius does
*not* transfer; the data-derived one does:

| Dataset          | median d1 | Otsu T | GMM T | auto-T (purity / recovery) | fixed T=80 recovery |
| ---------------- | --------- | ------ | ----- | -------------------------- | ------------------- |
| dino_dog_toy     | 71        | 83     | 72    | 0.97 / 0.94                | 0.94 *(80≈scale)*   |
| seoul_bull       | 131       | 112    | 103   | 0.995 / 0.98               | **0.67**            |
| seattle_backyard | 179       | 144    | 168   | 0.99 / 0.98                | **0.58**            |
| kerry_park       | 129       | 104    | 109   | 0.99 / 0.94                | **0.72**            |

The Otsu/GMM threshold lands in each dataset's own antimode and gives purity
0.97–0.998 with recovery 0.90–0.98 everywhere; the fixed `T=80` collapses off
dino. Neither estimator dominates — GMM tends a touch looser (higher recall),
Otsu a touch tighter (higher purity).

### Otsu's method

A classic image-binarisation threshold (Otsu, 1979). Given a 1-D histogram, it
chooses the split `t` that **maximises the between-class variance** of the two
groups (values ≤ `t` and > `t`) — equivalently, minimises the variance *within*
each group. For a split into classes with weights ω₀, ω₁ and means μ₀, μ₁,

```
σ_b²(t) = ω₀(t) · ω₁(t) · (μ₀(t) − μ₁(t))²
```

is maximised over `t`. Intuitively it finds the cut that makes the two sides as
internally tight and as far apart as possible — i.e. the **valley of a bimodal
distribution**. We run it on the `d1` histogram, where the two "classes" are
"has a real neighbour" (small `d1`) and "isolated" (large `d1`); the maximiser
is the antimode. It is cheap (one pass over the histogram) and uses no labels.

### Gaussian-mixture crossover

Fits a **2-component Gaussian mixture** to the `d1` values by EM
(`sklearn.mixture.GaussianMixture`), modelling the data as two overlapping
Gaussians — the "neighbour" mode and the "isolated" mode — each with its own
mean, variance, and weight. The threshold is the **crossover** between the two
means where the posterior responsibility flips from the low-mean component to
the high-mean one (a point becomes more likely "isolated" than "neighbour").

Unlike Otsu, which only sees the histogram and assumes a hard split, the GMM is
a soft probabilistic model that accounts for each mode's spread and size — which
is why it can land at a different (often slightly looser) point. We fit on
**log d1** (distances are positive and right-skewed, so the modes are more
Gaussian in log space) and exclude the `d1≈0` exact-duplicate spike, which would
otherwise capture a whole component.

## Files

- `sfm_descriptors.py` — descriptor-bank loader (above).
- `knn.py` — dependency-free chunked-GEMM exact k-NN (the oracle for validating
  the approximate index).
- `exp01_descriptor_kdtree_clusters.py` — clusters-vs-tracks probes.
- `exp02_ratio_threshold.py` — neighbour-distance structure + plots.
- `exp03_radius_clusters.py` — bounded-radius clustering, scored vs tracks.
- `exp04_auto_threshold.py` — data-derived radius (Otsu / GMM) + transfer check.
- `exp05_cluster_match.py` — the cluster matcher: writes a `.matches` file. Its
  canonical `--mode bgclusters` now calls the shipped production matcher
  (`sfmtool.background_floor_clusters` + `clusters_to_pair_matches`,
  `crates/sfmtool-core/src/cluster_match/`); the other exploratory modes still use
  the in-tree `sfmtool.KdForest` ANN index (`--exact` for the brute-force oracle).
- `exp07_cluster_vs_tracks.py` — scores clusters against the solve's tracks; its
  default `bgfloor` mode clusters with the production matcher, `global` keeps the
  POC numpy clustering as the comparison alternative.
- `exp06_solver_stability.py` — incremental/global solver seed sweep.
- `run_neighbor_pipeline.sh`, `run_both_solvers.sh` — end-to-end drivers
  (matches → COLMAP verify → solve → `sfm compare`).

The approximate NN index is now the in-tree randomized kd-tree forest
(`sfmtool.KdForest`, `crates/sfmtool-core/src/kdforest/`, spec
`specs/core/randomized-kdtree-forest.md`); the earlier standalone Rust
evaluation crate has been removed in favour of it.

Generated artefacts (`out/`) are gitignored; re-run to regenerate.
