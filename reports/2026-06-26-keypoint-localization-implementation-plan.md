# Implementation Plan: Keypoint-Localization Search Cache + Subpixel Refinement

_Date: 2026-06-26. An incremental plan for the two design specs that landed on
`main` in #135:_

- _**Grid search (Spec A):**
  `specs/core/keypoint-localization-search-cache.md` — the integer cross-view
  congealing search, accelerated by a render-once per-view cache, a hand-rolled
  AVX2 windowed-ZNCC kernel (centered-f32, then integer-`i16`), integer-tracked
  reads, and an `f32` search-resolution multiplier `m`._
- _**Subpixel refinement (Spec B):**
  `specs/core/keypoint-subpixel-refinement.md` — a standalone, continuous
  photometric (ECC / forward-additive Gauss–Newton) local refiner; the
  high-accuracy / ground-truth-approximating option._

_The prototype scripts on this branch (`scripts/exp_view_localization.py`,
`exp_normal_from_congealed_keypoints.py`, `exp_reference_refine.py`,
`_normal_strip_lib.py`) are the validated references these production
implementations should reproduce, and serve as oracles in the validation steps._

## Sequencing: A first, interleave a minimal B, gate the expensive work

**Recommendation: implement Spec A's *structure* first, interleave a *minimal*
Spec B early as the accuracy reference, then gate the costly optimizations of both
on a head-to-head measurement.** Not pure-sequential, not free-for-all interleave.

Why **A first**:

- **A is the producer.** B consumes A's output — the seed keypoint, the kept view
  set, and the consensus. There is little to build or test in B end-to-end until
  A's localization pipeline produces those.
- **A is partly landed.** The scalar correlation-accumulation `search_shift`
  (branch `optimize-search-shift`, ~2.5× kernel / ~1.9× embed, validated by
  `search_shift_matches_reference*`) is the foundation the cache and AVX2 build on.
- **A carries the larger, already-measured win.** On dino, localization is ~83% of
  `embed-patches`; within it `search_shift` ~56% and `render_context` ~33%. A's
  cache + AVX2 attacks both.

Why **interleave a minimal B early** (rather than finishing all of A first):

- The specs flag a real **architectural fork**: a *supersampled* grid (`m > 1`,
  Spec A) resolves sub-pixel directly and may be good enough, while the continuous
  LK refiner (Spec B) is the high-accuracy reference. Deciding between "grid-only"
  and "grid + LK" for the production path needs **both** in hand to measure — and
  the LK is the *reference* the grid's discrete accuracy is measured against.
- Answering that fork **before** sinking effort into the expensive optimizations
  (hand-rolled AVX2, `i16` stage 2, the LK accelerations) avoids over-investing in
  a path the measurement might not favor.

Why **not** pure interleave from step one: B's MVP still needs A's seed/consensus
to run on real data, so a minimal A precedes even the minimal B.

## Phases

### Phase 0 — Land the scalar foundation (done; merge)

The scalar correlation-accumulation `search_shift` on `optimize-search-shift`.
Already validated and rebased on `main`. **Action:** merge it. It makes the
integer search the accumulation form that the cache and AVX2 kernels extend.

### Phase 1 — Per-view cache + integer-tracked reads (Spec A, structural, scalar)

The render-once structural win, on the scalar kernel (no AVX2 yet).

1. **Render one expanded, frame-oriented cache per view** (planar, **centered**
   `f32`, plus a validity plane). Start with the unconditionally-correct
   `R + 4·search` size; add the `±search`-clamped `R + 2·search` variant behind
   validation.
2. **Round loop reads from the cache** — consensus cores *and* search candidates —
   dropping the per-round `render_context`. Split offset tracking into the integer
   read accumulator `iacc` and the parabolic **sub-pixel residual** (used only for
   convergence and as the hand-off seed; never fed back into reads, keeping every
   cache read integer-exact).
3. **`search_resolution_multiplier: f32`** on `KeypointLocalizeParams`, default
   `1.0` (the plan is the full-resolution grid; `m < 1` is a documented fallback).

- Files: `crates/sfmtool-core/src/patch/keypoint_localize.rs`;
  `keypoint_localize/prof.rs` (already there) for measurement.
- Validation: kept-view / registration **agreement** with the pre-change
  `search_shift_ref` congealing on the datasets (the documented sub-px behaviour
  change — *not* bit-equivalence); `embed-patches` point count stable; re-profile
  with `SFMTOOL_PROFILE=1` and confirm `render_context` collapses from
  per-(view, round) to per-view.
- Exit: integer localization fully on the cache; the render slice is gone; same
  kept-views/registrations within tolerance.

### Phase 2 — Minimal subpixel reference (Spec B, MVP — the accuracy reference)

The simplest correct continuous refiner — enough to *be* the high-accuracy
reference and to enable the measurement gate.

- New module `crates/sfmtool-core/src/patch/keypoint_subpixel.rs`:
  **forward-additive ECC Gauss–Newton**, per-view against a **single frozen
  (single-pass) consensus**, **bilinear** sampling from the **source pyramid**,
  gradient by the simplest correct route (the analytic bilinear gradient, or
  finite differences). Guard = accept only if the score improves, else keep the
  seed; out-of-frame / singular → keep seed. **Infinity patches refined too**
  (not skipped).
- PyO3 binding so the pipeline and tests can call it; optionally wire it into
  `_embed_patches.py` after `localize_keypoints`.
- Oracles: cross-check against `scripts/exp_reference_refine.py` /
  `exp_normal_from_congealed_keypoints.py` on this branch.
- Validation: synthetic recovery of a planted sub-pixel offset (< 0.02 px); patch
  **sharpness** + cross-view ZNCC on the datasets (sharpness should rise — the
  prototype's observation); the infinity case; guard correctness.
- Exit: a working, validated high-accuracy reference, callable from Python.

### Decision gate — measure the fork

With Phase 1 (grid, incl. the `m > 1` supersampled variant) and Phase 2 (the LK
reference) both available:

- Compare supersampled-grid sub-pixel (`m = 2, 3`) keypoints against the LK
  keypoints — the **accuracy gap** and the **speed** of each.
- Decide the production path: grid-only (`m > 1`) where its accuracy suffices,
  grid + LK where quality matters, or a selectable choice. **The LK stays the
  high-accuracy / ground-truth option regardless** — the gate only sizes *when*
  the cheaper grid is good enough.
- This decision directs where Phase 3 effort goes.

### Phase 3 — Optimize, guided by the gate

Build only the optimizations the gate justifies; validate and benchmark each
independently.

- **A-opt (if the grid is hot / on the path):** the hand-rolled **AVX2**
  centered-`f32` register-blocked kernel (with an `avx2_matches_scalar` test, like
  the fronto cache's `resample_avx2_matches_scalar`), then the **integer-`i16`**
  stage 2 (box-window normalization) gated on argmax agreement *and* a real
  speedup.
- **B-opt (if LK is on the production path):** the sampler value+gradient variants
  (`sample_bilinear_with_grad_u8`, `remap_aniso_with_grad`) + `WarpMap::get_jacobian`
  (the spec's "Design details" section); the incremental running-sum consensus with
  low-frequency IRLS weight refresh; inverse-compositional ECC only if its
  precompute amortizes against the refreshing consensus.

## Interleave vs. sequential — verdict

Mostly **sequential with one strategic interleave**. A's structural layer (Phases
0–1) goes first because it's the producer, it's partly done, and it's the bigger
win. A *minimal* B (Phase 2) is interleaved before the expensive optimization
work so the supersampled-grid-vs-LK fork can be measured and resolved at the gate.
The costly optimizations of both (Phase 3) are deferred until that measurement
says where they pay. Pure parallel interleaving from the start is not worthwhile —
B can't be exercised end-to-end without A.

## Risks / notes

- **Behaviour change, not equivalence.** The cache's integer-tracked reads build
  the consensus from integer-aligned cores; validate by kept-view/registration
  agreement, not bit-equivalence. (The search *kernel* itself stays equivalent —
  the existing `search_shift_matches_reference*` tests hold.)
- **Supersampled-grid cost grows ~`m²`** — the gate must weigh accuracy gained
  against that, not just whether sub-pixel is reachable.
- **Prototype parity.** The scripts on this branch are the validated references;
  reproduce their behaviour and reuse them as oracles rather than re-deriving.
- **Maturin rebuild.** Both specs touch `sfmtool-core` reached through PyO3;
  `pixi run maturin develop --release` before Python-side tests/profiling.
