# Grid Distortion Model (design exploration)

> _Status (2026-07-11): **idea-stage notes** — no production code, no
> production spec. We are exploring whether a grid-based distortion
> correction could support cluster-level geometric verification, using the
> small experiments in `scripts/exp_grid_distortion.py` to test the parts
> we are unsure about. Numbers below are from those experiments only and
> should be read as first evidence, not conclusions._

## The problem being explored

The setting: **a single physical camera capturing a scene from many
positions and angles** (e.g. a video sweep or photo set). Before any solve
exists we already have track clusters and, from cluster-patch refinement,
**affine correspondences** — per-member affine warps of photometric
quality. The question: **can the clusters with their affine
correspondences do a coarse bootstrap** — coarse cameras plus a coarse
camera correction — good enough to geometrically verify the clusters in
normalized ray coordinates (pixel-space constraints are wrong under
fisheye) and to seed a real solve, without depending on someone having
picked the right parametric distortion family?

The camera-correction piece being explored for that: a grid distortion
correction — a spline displacement field on a trivial base projection,
hopefully general enough to cover typical lenses without a model choice.
This is not a proposal to replace the parametric models elsewhere; COLMAP
interop, the solvers, and existing configs keep their models. (Captures
with several distinct cameras — fisheye video + high-res detail set +
drone — would carry one correction per intrinsics group, declared where
intrinsics already live: the per-directory `camera_config.json`; each grid
then pools observations from many frames, depths, and poses.)

The experiments below probe the camera-correction half of this
(representation size, joint observability with geometry, scheduling, real
tracks). **They do not yet touch the defining half — clusters and affine
warps as the bootstrap's input**; see the open questions.

## Model sketch

- **Base + delta.** The grid is a correction field on top of a trivial
  base projection — pinhole for narrow FOV, equidistant (`r = f·θ`) for
  fisheye — so "zero grid" is a sane camera and existing EXIF/config
  initialization is unchanged. Wide FOV works in the azimuthal angle chart
  `θ·(cos φ, sin φ)`, which stays smooth through and beyond 90°, where the
  pinhole chart diverges.
- **Cubic tensor-product B-spline** control grid over the observed image
  domain: C² smoothness, simple analytic Jacobians (fixed basis weights
  per observation), cheap to evaluate.
- **Vector-field symmetry is covariant**: under a mirror in x, the x
  component is odd and the y component even (`d(−x, y) = (−dx, dy)`),
  about the principal point (which stays an explicit K parameter).

## Coarse → fine, along several dimensions

The hope is that typical lenses resolve at the coarse levels and each
refinement has to be earned by the data.

1. **Symmetry class**: `L0` radial 1-D profile `g(r)·r̂` → `L1`
   quadrant-symmetric 2-D grid (4-fold reflection tying, ¼ the DOF; also
   removes the skew and principal-point-shift gauge modes by construction)
   → `L2` full 2-D grid (decentering/tilt terms need this level).
2. **Control-grid resolution**: knot doubling within a level
   (2→4→8→16 per axis), initialized by refinement of the coarser fit.
3. **Solve curriculum**: base K first, then levels/resolutions unlock in
   order, with smoothness regularization annealed per level.
4. **Radial extent, center-out**: solve with central observations first —
   where the candidate models (including "none") agree and distortion is
   weakest — then admit growing radii, unlocking `L0` knots outward in
   step. Plausible benefits: less sensitivity to the initial model guess,
   data admission ordered by measurement reliability (peripheral
   detections are blurrier and more foreshortened), and it composes with
   the knot-resolution dimension. Experiment C suggests what it does *not*
   do: it did not help against absorption (see below).
5. **Feature size, large → small**: truncating to the largest features is
   approximately extracting from a down-rezzed image, so it is a data-side
   resolution pyramid that costs nothing — `.sift` files are already
   sorted by descending feature size, so "top-K rows" is the whole
   implementation. Coarse stages can run on large features only
   (fewer, more stable, cheaper), admitting smaller features as the model
   refines. (Noted as a dimension; not yet exercised by the experiments.)
6. **Per-group unlocks, gated on coverage *and* geometry strength**: a
   level/region unlocks for a group only when its grid cells have enough
   observations *and* the group's view diversity makes distortion
   observable at all — in experiment B a frontal-only capture absorbed
   the field into bent structure regardless of coverage. Unlock criteria
   should be data-derived (held-out residual improvement, ray-diversity
   per cell), not fixed constants.

**Gauge constraints** (likely needed at every level): zero displacement
and zero linear term at the center — a linear field is indistinguishable
from changing `f`/PP/skew, and the isotropic-scale mode survives even at
`L0`. Within a group, the grid would hold the (stable) lens geometry while
per-capture focal breathing stays in `f`.

## Experiments (`scripts/exp_grid_distortion.py`)

- **A — expressiveness.** For known parametric models (the Kerry Park
  OPENCV_FISHEYE, a typical SIMPLE_RADIAL, an OPENCV Brown–Conrady with
  tangential terms): how many DOF per level until the parametric field is
  matched to well under solver noise (~0.05 px)? Does the antisymmetric
  part force `L2`, and how large is it?
- **B — absorption/degeneracy.** Synthetic BA (poses + points + `f` + full
  grid, known ground truth): does the grid absorb structure or `f`?
  Conditions: frontal vs orbit geometry, full-field vs center-only
  coverage, gauge constraint on/off.
- **C — center-out curriculum.** Does staged radial admission help the
  weak-geometry case?
- **D — real data.** COLMAP South Building: re-solve a subset of the real
  observations with pinhole + grid and compare against the dataset's
  solved SIMPLE_RADIAL model.
- **E — the bootstrap itself** (`scripts/exp_cluster_bootstrap.py`): from
  our clusters and cluster-patch warps alone — no solver initialization —
  alternate a missing-data affine factorization (per-image `[M_i | t_i]`,
  per-cluster `X_c`, ALS) with an `L0` radial-profile fit, on a window of
  a real capture; score against a reference solve of the same subset.

### Results (2026-07-11 run)

**A — expressiveness.** DOF = free parameters after symmetry tying; errors
are px over the sampled field (fisheye in angle-chart px).

- **Kerry Park OPENCV_FISHEYE** (field rms 9.4 px, max 14.4 px): the
  radial profile fits efficiently — `L0` reaches 0.14 px rms with 8 knots
  and 0.05 px rms / 0.40 px max with 12 knots, while a 2-D grid needs
  288 DOF to reach 0.11 px rms. `L1` and `L2` produce identical fits here
  (the field is fully symmetric — a sanity check).
- **SIMPLE_RADIAL** (k = −0.12): exact at `L0` with the minimum 4 knots
  (the displacement is cubic in `r`, which a cubic spline represents
  exactly — another sanity check more than a result).
- **Brown–Conrady with tangential terms** (k1 −0.28, k2 0.07, p1 0.0012,
  p2 −0.0007; field rms 26.8 px): `L0` and `L1` plateau at the
  antisymmetric floor (0.4897 px rms, max 1.62 px — the decentering
  component, which no symmetric level can represent), and `L2` clears it:
  0.044 px rms at 72 DOF, 0.009 px at 128 DOF.

In these fits, radial-first ordering looks right, the antisymmetric floor
matches the `L1`→`L2` transition, and the parametric models tested are
matched with ≤ 128 free parameters.

**B — absorption/degeneracy.** Synthetic BA, 20 cameras × 400 points,
0.3 px noise, GT SIMPLE_RADIAL k = −0.10, full 6×6 grid + f + poses +
points all free, f started 2% off. "Map-shape err" compares the learned
composite projection `f·x + grid(f·x)` to the GT map at observed
locations, **modulo one uniform scale** — two gauges have to be modded
out: the linear-gauge split (with the grid's linear mode pinned, `f`
carries the distortion's best-linear part) and the map-scale mode (an `f`
change compensated by cameras retreating / a depth rescale; it appears to
be the weakest-observable distortion mode in these runs, so shape recovery
is judged without it, keeping in mind that the mode itself is only weakly
determined).

| geometry | coverage | gauge | reproj px | struct rmse | map-shape err rms / max px |
|----------|----------|-------|-----------|-------------|-----------------------------|
| frontal  | full     | on    | 0.298     | 0.088       | 1.39 / 9.49                 |
| frontal  | full     | off   | 0.298     | 0.082       | 1.30 / 9.43                 |
| frontal  | center   | on    | 0.268     | 0.118       | 0.18 / 0.78                 |
| orbit    | full     | on    | 0.288     | 0.0047      | 0.27 / 1.48                 |
| orbit    | full     | off   | 0.290     | 0.0029      | 0.15 / 0.78                 |
| orbit    | center   | on    | 0.282     | 0.0115      | 0.19 / 0.50                 |

In this setup the deciding factor was geometric diversity, not coverage
and not the gauge: the frontal scene (slab + translating cameras) bends
the structure and loses much of the field's higher-order shape — the
classic dome degeneracy, which the gauge constraint does not address (it
fixes parameter interpretation, not the radial↔depth trade). The orbit
scene keeps the structure nearly intact and recovers the map shape to
0.15–0.27 px rms, even when only the sensor center is covered. Whether
these synthetic scenes bracket real captures well is an open question.

**C — center-out curriculum vs joint solve** (frontal scene, gauge on):
staged solving — central observations with the grid frozen (r < 0.35),
then r < 0.65 with the grid unlocked, then everything — tested as a
possible absorption defense. In this experiment it was not one:

| solve | struct rmse | map-shape err rms / max px |
|-------|-------------|-----------------------------|
| joint (one shot)          | 0.073 | 1.34 / 7.20 |
| stage 1 (r<0.35, no grid) | 0.057 | 3.42 / 18.9 (grid still zero) |
| center-out final          | 0.069 | 1.30 / 7.78 |

A plausible reading: the dome ambiguity is scale-continuous — plain BA on
the central data alone already bends the structure by a smaller version of
the same trade, and later stages warm-start inside the same valley, ending
where the joint solve ends. Center-out may still be worth having for its
other properties (dimension 4 above), but on this evidence it is not an
absorption defense.

**D — real data (COLMAP South Building).** 128 images around a building,
one shared solved SIMPLE_RADIAL camera (f 2559.7, k = −0.0205, 3072×2304).
A 24-image / 800-point / 3406-obs subset of the real keypoints, re-solved
from a pinhole start (f off by 2%, grid zero, poses/points from the
reconstruction):

- reference model's own reprojection rms on the subset: 0.770 px; the
  grid solve reaches 0.450 px. The grid is fitting residual structure the
  one-parameter model doesn't — though with 72 grid DOF vs 1, some of
  that could be overfitting; a held-out protocol would tell.
- camera-map shape vs the solved model: 1.27 px rms / 9.7 max over a
  4.1 px rms (18.4 max) reference field, worst at the sparsely-observed
  footprint corners; structure drift vs the reference reconstruction
  0.0075 (similarity-aligned). Some of the 1.27 px may be the
  SIMPLE_RADIAL reference's own error rather than the grid's — this
  comparison cannot distinguish them.
- `f` stayed at its perturbed init (+2%), compensated by geometry, and
  the fitted map scale (0.976) brought the maps into agreement —
  consistent with the map-scale mode being soft on real data too.

**E — coarse bootstrap from real clusters (first attempt).** DinoLedge
window: 40 frames at stride 3 (~4 s of a walking capture), full cluster
pipeline (sift → `match --cluster` → `cluster-patches`), plus an
incremental solve of the same subset as reference (SIMPLE_RADIAL,
f 2749, k1 +0.0135 — a nearly distortion-free, likely in-camera-corrected
phone lens). 3374 clusters span ≥ 6 of the 40 frames. Bootstrap =
positions-only affine factorization alternating with an `L0` fit, from
scratch; camera comparison is gauge-aligned (the factorization's global
3×3 affine ambiguity is solved for before measuring row-space angles).

| clusters | obs | factorization residual rms | camera row-space angle mean |
|----------|------|---------------------------|------------------------------|
| 150      | 1082 | 2.12 px | 16.5° |
| 500      | 3604 | 2.84 px | 19.5° |
| 3374     | 24145 | 3.22 px | 9.6° |

- **Coarse cameras from clusters alone appear feasible**: ~10° mean
  attitude error with the full cluster set, no solver anywhere. Whether
  10° is good enough to verify clusters and seed a solve is the next
  question, not answered here.
- **Distortion recovery is untestable on this dataset**: the reference
  profile over the window footprint is only 0.07 px rms (mod linear) —
  far below the 2–3 px factorization residual. The bootstrap "recovered"
  a 0.6–0.7 px profile, i.e. spurious absorption of affine-model error
  into the radial term. Lesson: the profile fit needs a significance gate
  against the factorization noise floor, and testing the distortion half
  of the bootstrap needs a genuinely distorted capture (South Building
  imagery through our cluster pipeline, or the Kerry Park fisheyes).
- **The warp check was inconclusive for the same reason**: the
  multi-cluster warp rank-2 misfit did not change between the zero and
  estimated profiles (~0.13 both) — there was no distortion signal to
  sense. What the affine correspondences buy remains untested.
- A negative control worth keeping: sampling the *whole orbit* at stride
  24 instead of a window broke the factorization outright (~60° camera
  angles even after gauge alignment) — the affine model needs windows of
  nearby views, as the design assumed; chaining windows is the untested
  next mechanism.

**E across datasets (same run date).** Four captures through the same
pipeline (window subsets, full cluster pipeline, reference solve where
possible). dino / seoul / seattle share one physical phone (dino = photo
mode, the others = video mode); Kerry Park is one fisheye of the rig pair,
scored against the rig calibration instead of a solve, with the
observation footprint cut at 0.45·rmax (the perspective chart diverges
toward a fisheye's edge) and the profile seeded from the config focal
(the production base-plus-config-init setting).

| dataset | images | clusters (max run) | resid rms | cam angle | ref profile scale | profile err (mod linear) |
|---------|--------|--------------------|-----------|-----------|-------------------|--------------------------|
| dino (photos)   | 15 | 2916 | 0.98 px | **2.7°** | 3.09 px | 2.85 px (absorbed) |
| seattle (video) | 26 | 611  | 0.65 px | **3.3°** | 11.22 px | 11.75 px (absorbed) |
| seoul (video)   | 17 | 26   | 0.27 px | **6.4°** | 0.08 px | — (no signal) |
| kerry (fisheye, seeded) | 24 | 446 | 2.48 px | n/a | 16.73 px | **5.65 px** |
| kerry (no seed, ablation) | 24 | 446 | 0.79 px | n/a | 16.73 px | 16.64 px (fully absorbed) |

Readings:

- **Coarse cameras are consistently cheap**: 2.5–6.5° mean gauge-aligned
  attitude error across every phone dataset, down to 26 clusters on tiny
  270×480 seoul frames.
- **Positions-only distortion recovery fails on window-sized arcs, on
  real data**: dino's 3 px and seattle's 11 px profiles were absorbed
  into cameras/structure at sub-pixel residuals, and kerry's unseeded
  ablation absorbed a 16.7 px fisheye deviation *while improving* the
  fit residual — experiment B's absorption, reproduced on real captures.
  The config-seeded kerry run kept roughly two thirds of the profile,
  so the base + config-init design earns its place, but window geometry
  alone cannot pin the profile.
- **The warp misfit metric needs work before it can arbitrate**: on
  kerry it dropped sharply with the estimated profile (0.134 → 0.046–
  0.065) in *both* the seeded and fully-absorbed runs, so it is
  responding to something other than profile correctness (possibly any
  radially-contracting Jacobian); treat the "warps sense distortion"
  signal as suggestive, not established.
- Same-phone footnote: the captures from one phone show wildly different
  effective distortion — photo mode (dino, ~3 px profile), full-res video
  (DinoLedge, ~0.07 px), downscaled test videos (seattle ~11 px, seoul
  ~0.08 px). The phone has multiple physical cameras (wide-angle vs
  normal), so this is likely a mix of lens selection (seattle's strong
  barrel reads as the wide-angle) and per-mode processing/corrections —
  we cannot separate the two from this data. Either way the practical
  rule holds: corrections attach to the capture group, never the device.

**F — can a pinhole start evolve into a fisheye, center-out?**
(`scripts/exp_pinhole_to_fisheye.py`, Kerry Park left fisheye, 24 frames,
446 clusters.) Observations admitted in growing radial rings; per ring,
plain pinhole factorization (control) vs factorization alternating with a
warm-started `L0` profile refit (the center-out evolution); at the end, a
`s·f·tan(r/f)` fit to the evolved profile attempts the model-class
upgrade. Reference-field scale grows violently with the ring (mod-linear
rms 2 px at ρ=0.25 → 60 px at 0.55 → 269 px at 0.65, past θ=90°).

| ring ρ | ref field rms | pinhole resid | evolve resid | profile recovered |
|--------|---------------|---------------|--------------|--------------------|
| 0.25 | 2.0 px  | 1.05 px | 1.00 px | ~0 |
| 0.45 | 19.0 px | 0.73 px | 0.78 px | ~0 |
| 0.55 | 60.4 px | 1.15 px | 1.18 px | ~1 px of 60 |
| 0.65 (θ>90°) | 269 px | 0.92 px | 1.03 px | ~0 |

The answer for **this mechanism** is no. Two qualifications matter:

- **Trimming confound / contamination floor.** The table's residuals are
  over the robust fit's surviving ~55% of observations. Re-running with
  trimming off: even at the *central* ring, where the model class is
  unquestionably adequate, the median residual is 6.3 px and only 28%
  of observations sit within 3 px — about half the cluster data is
  contamination/noise independent of any model question. Against that
  floor, growing the ring to a 269 px field moved the no-trim median
  only to ~8 px: the fisheye contributes almost nothing to the residual
  distribution at any radius.
- **Why the absorption happens** (the mechanism, not magic): (1) the
  factorization sees only cross-frame consistency, never absolute
  positions — lens distortion is a *static* warp, and each cluster's
  free 3D point absorbs the static component exactly (with zero camera
  motion, any distortion would be perfectly invisible); (2) distortion
  is therefore visible only through *motion* — the warp's local Jacobian
  modulates apparent motion, giving a signal of order (Jacobian
  variation) × (within-window motion) ≈ a few px here, not the 60–269 px
  field; (3) the affine model's nuisance freedoms span exactly that
  signal — per-cluster depth rescales parallax, the free 3D direction
  tilts implied motion (an affine model has no epipolar constraint), and
  the per-image affine eats anything global. What remains is below the
  contamination floor.

**Scope of the negative:** F rules out position-only, affine, single
window. It does *not* rule out within-window lens signals via: the
**warps** (an affine correspondence measures the map's Jacobian between
two image locations — the static-absorption argument structurally does
not apply to it; a central-reference/peripheral-member Kerry cluster
carries the radial compression in the warp itself, motion or no motion);
**epipolar rigidity** (perspective two-view constraints are stiffer than
affine — fisheye periphery violates pinhole epipolar geometry
measurably); **wider baselines**; and **cross-window pooling** (bent
per-window worlds cannot agree on one shared camera). Center-out radial
admission survives as an *admission schedule* (consistent with
experiment C), not as an information source.

**Working conclusions (to be revisited):**

0. **The coarse bootstrap should carry no distortion DOF at all** (user
   direction after the cross-dataset runs, and the data supports it):
   at window scale, distortion is indistinguishable from — and absorbed
   by — camera/structure DOF, with no measured harm to the coarse
   cameras (2.5–6.5° everywhere, including a fully-absorbed 16.7 px
   fisheye deviation at the *best* residual of the set). Estimating it
   early is therefore both unnecessary and unidentifiable (it only
   manufactures overfit, e.g. DinoLedge's spurious 0.6 px profile).
   The base *chart* from config stays (raw fisheye pixels are not even
   approximately affine beyond ~60° — that's domain validity, not a
   parameter), and grid DOF enter only later, once pooled geometry
   passes the strength gate — which a single window never would.
   Corollary for the verifier: with pinhole-only coarse geometry, inlier
   residuals carry a systematic radius-dependent component, so the
   inlier/outlier split must be data-derived per window, not a fixed
   pixel gate.

1. The coarse-to-fine curriculum probably needs a geometry-strength gate,
   not just a coverage gate: distortion DOF should unlock only when the
   group's view diversity makes them observable (possible proxies:
   angular/baseline diversity of rays per grid region, or the Fisher
   information of the radial-depth mode).
2. Coverage gating still matters for *where* the grid is trusted — the
   recovered map is only meaningful over the observed footprint.
3. The gauge constraint seems worth keeping for identifiability and
   stable `f` semantics (f ≈ best-linear scale over the coverage), not as
   an absorption defense.

## Sketch: how this could sit in the code

Current thinking, if the idea survives further scrutiny: **two new
variants on the existing camera-model enum**
(`sfmtool_core::camera::intrinsics::CameraModel`, whose variants mirror
the COLMAP models and serialize through the same `SfmrCamera` parameter
convention), split by the base projection's chart:

- **`PinholeGrid`** — perspective chart, base pinhole. Zero grid ≡
  `Pinhole`. Narrow/medium FOV.
- **`FisheyeGrid`** — azimuthal-equidistant chart (`r = f·θ`), base
  equidistant. Zero grid ≡ ideal equidistant fisheye; usable through and
  beyond 90°.

As enum variants they would be ordinary cameras wherever sfmtool-native
code consumes a `CameraModel` — remap/undistort, frustums, the patch
pipeline, `.sfmr` files, and the verification path that motivates them.

No "base radial" variant is planned: experiment A suggests the `L0` radial
profile covers realistic radial models with ≤12 knots, so a polynomial
base under the grid would double-parameterize the same function space. A
known parametric model would instead be an *initialization source*: fit
the grid to its field once and bake it into the control values.

Per model, the stored parameters would be
`[fx, fy, cx, cy] + hierarchy descriptor (level, knots per axis) +
control values` — variable-length params vectors, which the
`camera_config.json` / `.sfmr` / `.camrig` model+params pattern already
accommodates. `L0` at minimum knots with zero values is the base model, so
the hierarchy nests in one name and the stored form doubles as the
curriculum state.

**Conversions at the boundary.** pycolmap and the COLMAP DB never speak
grid, so both directions need supported conversions:

- *Into the grid*: a known parametric model (config, EXIF inference, or a
  prior solve) seeds the grid by fit-and-bake.
- *Out of the grid*: once a grid is solved, export paths convert to the
  parametric variant that fits best over the observed footprint (recording
  the fit residual), or undistort to `Pinhole` — e.g. `to-colmap-db` and
  any solver hand-off would pick up the converted model.

The parametric variants remain the interop and solver formats.

**Direction:** probably store the *projection* correction (ideal→observed
— what reprojection evaluates constantly) and unproject with a few Newton
steps on the C² field. Revisit if unprojection shows up hot in a profile.

## Open questions

- **The bootstrap's distortion half.** Experiment E showed coarse cameras
  from clusters alone, but its dataset had no distortion to recover; the
  same experiment needs a distorted capture (South Building images run
  through our cluster pipeline, or the Kerry Park fisheyes with a
  `FisheyeGrid`-chart variant), plus a significance gate so the profile
  fit cannot absorb affine-model error.
- **What the warps add.** Still untested (E's warp check had no signal to
  sense). Point positions constrain the distortion field; an observed
  affine warp involves the distortion map's *derivative* at both image
  locations. Whether those pointwise gradient constraints help — and
  whether they weaken the dome degeneracy of experiment B — is unknown.
- **Window chaining.** E worked on one window; a whole capture needs
  windows chained into a consistent coarse geometry (and the full-orbit
  negative control shows windows are mandatory, not optional).
- **Is ~10° coarse-camera accuracy sufficient** to verify clusters and
  seed a solve? Needs the verifier prototype to answer.
- Whether the azimuthal chart suffices for the widest lenses we care
  about, or a direction-sphere parameterization is eventually needed.
- A held-out–reprojection evaluation protocol (experiment D shows
  model-vs-model comparison confounds reference error with grid error,
  and can't rule out overfitting).
- Whether the synthetic frontal/orbit scenes bracket real capture
  geometries well enough to trust the absorption conclusions.
