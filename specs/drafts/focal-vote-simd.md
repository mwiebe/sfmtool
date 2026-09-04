# Draft: SIMD pass over the focal-vote kernel's residual loops

Change proposal. Converts to nothing on its own -- the result of this work is
faster code with bit-identical output, plus updated implementation notes in
[`focal-vote.md`](../core/geometry/focal-vote.md) if any code structure worth
documenting changes.

## Goal

Cut the per-core cost of the structure-free focal vote's hot loops with
lane-parallel SIMD, under a **bit-identity bar**: `FocalVoteResult` (and
therefore every `estimate_intrinsics` product and every seed annotation) must
be bit-for-bit unchanged, SIMD on or off, on every input.

## Measured profile (what to target, and what not to)

Single-thread (`RAYON_NUM_THREADS=1`), release build, phase accumulators via
`Instant` around each function, minimum-of-warm-run. Two captures:
`kerry_park` (39,481 obs, 480x480, fisheye) and `dino_dog_toy` (409,838 obs,
2040x1536, pinhole).

Pinhole-only vote (the Auto policy's screening pass):

| phase | kerry | dino | calls (kerry) |
|---|---|---|---|
| wall | 709 ms | 150 ms | |
| `homography_dlt` | 608 ms | 102 ms | 249,444 |
| `score_h` | 78 ms | 6 ms | 249,510 |
| `epipolar_pair_outcome` (whole F side) | 19 ms | 21 ms | 18 |

Both-columns vote (an escalated / `--model fisheye` run) adds the scan cells;
their split is stable across the two captures:

| phase | kerry | dino | calls (kerry) |
|---|---|---|---|
| wall | 2952 ms | 2778 ms | |
| `null_from_rows` | 1384 ms | 1350 ms | 540,426 |
| `rotation_residuals` | 327 ms | 488 ms | 254,082 |
| `epipolar_residuals` | 261 ms | 456 ms | 540,426 |
| `kabsch` | 114 ms | 156 ms | 254,082 |
| ray rebuild loops | 25 ms | 38 ms | 8,742 |
| `epipolar_rows` | 2 ms | 4 ms | 4,130 |

Two conclusions frame the scope:

1. **The kernel is eigen-bound, not residual-bound.** `null_from_rows` and
   `homography_dlt` spend their time in nalgebra's 9x9 `symmetric_eigen`, run
   once per minimal-sample hypothesis (the RANSAC/scan inner step). That is
   ~50% of an escalated run and ~85% of kerry's pinhole-only run, and SIMD
   cannot touch it under the bit-identity bar. It is named as a follow-up
   below, out of scope here.
2. **The addressable surface is the residual loops** --
   `epipolar_residuals`, `rotation_residuals`, `score_h` -- roughly 0.9 s of
   an escalated run's 2.8-3.0 s serial. The realistic whole-kernel win for
   this pass is therefore on the order of 15-25% on escalated runs; report
   the measured number, do not promise more. The batched-trig idea for the
   ray rebuilds is refuted by measurement (25-38 ms) -- skip it.

## The bit-identity rule

A lane-parallel evaluation is bit-identical exactly when each lane performs
the same IEEE-754 operations in the same order as the scalar code did for
that element. Concretely:

- **Vectorize across points** (lane = one correspondence), never across the
  arithmetic of one point. Elementwise add / mul / div / sqrt / abs / min /
  max / compare are IEEE-exact per lane and thus safe.
- **No vectorized horizontal reductions.** A sum whose order changes is a
  different result. Per-point dot products of fixed length 3 keep their fixed
  scalar order inside each lane (mul, mul, mul, add, add -- same sequence).
- **No vectorized transcendentals.** `acos` / `sin` / `cos` are libm calls; a
  polynomial SIMD replacement is not bit-identical. They stay scalar.
- Where a guard produces `INFINITY` (the `symmetric_transfer_sq` degenerate
  branch), the SIMD path reproduces it with compare + blend, not by skipping
  the lane.

`kabsch`'s covariance accumulation and `homography_dlt`'s AtA accumulation
are order-sensitive sums over points -- out of scope (their cost is mostly
the eigen/SVD anyway).

## Targets, in order

1. **`column_scan::epipolar_residuals`** (~260-460 ms): per point,
   `n = E*r1[i]` (or transposed), `out[i] = min(|n . r2[i]| / max(|n|, 1e-15), 1)`.
   All exact ops, ~25 flops/point, 4-wide f64 AVX2. The rays live in
   `Vec<Vector3<f64>>` (24-byte interleaved); either pack lanes with
   scattered loads or restructure to SoA -- implementer's choice, measured.
   If SoA wins, the arrays are built once per grid point in `scan_epipolar`
   and shared with `fit_epipolar`; keep `epipolar_rows` reading whichever
   layout is chosen without a second copy of the data.
2. **`column_scan::rotation_residuals`** (~330-490 ms): per point,
   `(rot * r1[i]) . r2[i]`, clamp, `acos`. Split: vectorized matvec + dot +
   clamp into a scratch buffer, then a scalar `acos` tail loop. Measure the
   acos share first (add a temporary sub-timer); if acos dominates the loop,
   say so in the report and take only the matvec win.
3. **`homography_estimation::score_h`** (~6-80 ms here, but a shared
   primitive -- every `estimate_homography` caller in the workspace
   benefits): `symmetric_transfer_sq` over all points against one `H`,
   `H_inv`. All exact ops including the degenerate-guard blend.

Each target keeps the scalar implementation as the dispatch fallback (it IS
the current code, moved verbatim), with the kdforest idiom:
runtime `is_x86_feature_detected!("avx2")`, `#[target_feature(enable = "avx2")]`
`unsafe` inner kernels, and a kill switch
`SFMTOOL_FOCAL_VOTE_NO_SIMD` (mirroring `SFMTOOL_KDFOREST_NO_SIMD`,
`distance.rs:186`). SSE2 tiers are optional -- add only where the AVX2 win
justifies the second kernel.

## Verification

- **Bit-identity tests** per target, kdforest-style
  (`kdforest/distance/tests.rs`): scalar vs SIMD on seeded random inputs
  (including lengths around the lane width and the tail) compared with
  `f64::to_bits`, plus the degenerate-guard case for `score_h`.
- The existing focal-vote and estimate-intrinsics test suites must pass
  unchanged -- they pin the kernel's outputs and stay green by construction
  if the rule above is followed.
- **End-to-end bit-identity** on the two observation dumps: full
  `FocalVoteResult` (all floats hex-compared) for pinhole-only and
  both-columns, SIMD on vs `SFMTOOL_FOCAL_VOTE_NO_SIMD=1`.
- `pixi run cargo fmt && pixi run cargo clippy --workspace && pixi run doc`;
  `pixi run cargo test --workspace`. No binding surface changes, so no
  maturin/pytest leg is required unless the agent touches `sfmtool-py`.

## Bench protocol

Observation dumps live at the session scratchpad under `simd-profile/`
(`kerry_park.bin`, `dino_dog_toy.bin`; format: `u64 n, u32 w, u32 h,
u32 ci[n], u32 ii[n], f64 xy[2n]`, little-endian -- regenerate from the
scratch workspaces' `clusters.matches` via
`sfmtool._commands.estimate_intrinsics._load_observations` if needed).
Drive with a scratch `cargo` example (delete before finishing) that loads a
dump and times `focal_vote_with_options` for `[Pinhole]` and
`[Pinhole, EquidistantFisheye]`: warm pass, then min over 3 timed passes.
Report per-target and whole-kernel numbers at `RAYON_NUM_THREADS=1` and at
full threads, SIMD on vs off, both captures.

## Out of scope, named for the follow-up decision

The eigen-bound minimal solvers are the real ceiling: a fixed-size
elimination null-space (8x9 for `EPI_SAMPLE = 8`, and a dedicated 4-point
homography solve replacing the per-iteration normalized DLT + 9x9 eigen)
would attack ~50-85% of the kernel where SIMD attacks ~30%. That change is
**output-changing in the last bits** (the null vector of a noisy minimal
sample differs in fp noise between eigen and elimination), so it trades the
bit-identity bar for the fleet A/B protocol: verdicts, confirmations and
consensus focals invariant across all 42 fleet entries, spread/mass drifts
reported. Do not start it under this draft; it is a separate decision.
