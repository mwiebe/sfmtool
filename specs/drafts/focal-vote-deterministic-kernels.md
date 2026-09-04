# Draft: deterministic fast kernels for the focal vote's minimal solvers and rotation residuals

**Status:** Draft

Change proposal. Two hot paths of the structure-free focal vote
([`focal-vote.md`](../core/geometry/focal-vote.md)) are replaced by
self-contained deterministic arithmetic: the rotation-residual `acos` becomes a
fixed polynomial evaluated identically by the scalar and SIMD paths, and the
minimal-sample null spaces (the essential cell's 8-row samples, the
homography's 4-point DLT) become a direct 8×9 elimination instead of
`AᵀA` + 9×9 `symmetric_eigen`. Measured on the two profiling captures:
pinhole-only votes 2.1–4.9×, escalated both-columns votes 2.2–3.1× serial and
1.7–2.4× at full threads.

## Why the invariant moves, and to what

Both changes alter outputs in the last bits, so the byte-parity bar of the
SIMD pass does not apply. What replaces it is stronger where it matters:

- **Scalar ≡ SIMD by construction.** The polynomial is exact IEEE arithmetic;
  its scalar twin performs the same operations in the same order, so the two
  dispatch arms produce identical bits and `SFMTOOL_FOCAL_VOTE_NO_SIMD`
  remains a pure performance switch.
- **Platform determinism.** `f64::acos` calls the platform libm, so today's
  vote already differs in last bits across operating systems. The polynomial
  and the elimination are ordinary `f64` arithmetic — the changed paths
  compute the same bits on every platform and architecture. The scalar
  polynomial is therefore **not** gated on `x86_64`: every platform uses it.
- **Acceptance is fleet verdict-invariance**, run by the operator on the
  42-capture fleet: every camera-model verdict, confirmation and escalation
  list identical between the old and new kernels, consensus-focal drift
  reported. The experiment this draft productizes measured zero hard
  differences, worst consensus-focal relative drift 2.4e-14, worst
  diagnostic-leaf drift 1.5e-9.

Two diagnostic environment flags restore the old arithmetic for forensics
(they are not compatibility switches and appear in no production path):
`SFMTOOL_FOCAL_VOTE_LIBM_ACOS` restores the libm `acos`,
`SFMTOOL_FOCAL_VOTE_EIGEN_MINSOLVE` restores the eigen minimal solvers.

## Polynomial `acos`

`acos(d)` for `d ∈ [-1, 1]` via the asin core, branch-free in the SIMD form:

```text
a = |d|;  big = a > 0.5
z = big ? (1 − a)/2 : a²          (z ∈ [0, 0.25])
s = big ? √z        : a           (s ∈ [0, 0.5])
asin(s) = s + s·z·P(z)            (P = degree-13 polynomial, Horner from A13)
acos(d) = big ? (d < 0 ? π − 2·asin : 2·asin)
              : π/2 − copysign(asin, d)
```

The one coefficient table is shared by the scalar and vector evaluations (a
second copy would be free to drift — the reason `numeric.rs` exists):

```text
A0  =  0.16666666666666666     A7  =  0.01158174875867126
A1  =  0.07500000000000406     A8  =  0.009513962068006003
A2  =  0.04464285714150171     A9  =  0.009846417300000855
A3  =  0.030381944571586314    A10 =  0.0012909120006875199
A4  =  0.02237215339997836     A11 =  0.02336097864008306
A5  =  0.017352913993464503    A12 = -0.024103731139602444
A6  =  0.013962288953824856    A13 =  0.03238761648605816
```

Derivation: degree-13 Chebyshev interpolation of
`P(z) = (asin(√z) − √z)/(z·√z)` on `[0, 0.25]`, node values computed with
exact rational series arithmetic, converted to the monomial basis. Measured
accuracy: max 1 ULP against libm over dense `[-1, 1]` sampling plus
adversarial near-`±1` populations; the trailing coefficients are
interpolation artifacts that compensate one another and are exact as written.

Placement: the coefficient table and the scalar evaluation live in
`geometry::numeric` (platform-independent); `geometry::simd::acos_pd` is the
AVX2 form reading the same table, with blends in place of the branches
(`copysign` as `or(r, and(signmask, d))`; the big/small and sign selections as
compare + `blendv`). `rotation_residuals` uses the polynomial in **both**
dispatch arms — the scalar fallback, the vectorized body, and the ragged tail.

NaN propagates: a NaN cosine (already possible only from NaN rays) yields a
NaN residual through either form, as libm did.

## Elimination minimal solvers

`numeric::null9_from_8rows` computes a unit-norm right null vector of an 8×9
system by Gaussian elimination with partial pivoting (strict `>` comparison,
first maximal pivot kept — the tie rule is part of the determinism contract).
A generic rank-8 sample has a one-dimensional null space; elimination returns
it directly, avoiding both the `AᵀA` squaring (which squares the condition
number) and the iterative eigen decomposition. Rank-deficient input takes the
last free column with the other free coordinates zero — a deterministic
member of the null space; such degenerate samples score few inliers and lose
the RANSAC regardless. A design with no pivot at all, or a non-finite result,
is `None`.

Call sites, all minimal-sample only:

- `fit_epipolar`'s hypothesis loop (samples of exactly `EPI_SAMPLE = 8` rows)
  through a thin wrapper that reshapes the 9-vector row-major into the `3×3`
  epipolar matrix, replacing `null_from_rows` there.
- `homography_dlt` when `N = 4`: the eight DLT rows of the Hartley-normalized
  points feed the same primitive; denormalization and unit-Frobenius scaling
  are unchanged. (`E` and `H` consumers are scale- and sign-invariant, so the
  normalization convention is cosmetic.)

**Unchanged**: the consensus refits (`null_from_rows` over an arbitrary index
set is a least-squares smallest-direction problem, not an exact null space),
`kabsch`, `estimate_fundamental`, and `relative_pose` (its own consumers
carry their own validation; widening the scope there is a separate decision).

## Tests

- `acos` accuracy: absolute error against `f64::acos` bounded by `5e-16`
  over a dense `[-1, 1]` grid plus near-`±1` samples; exact behavior at
  `±1`, `±0.5`, `±0.0`; NaN in, NaN out.
- Scalar/SIMD bit-identity: `acos_poly_scalar` vs `acos_pd` compared with
  `to_bits` over seeded inputs including the boundary and NaN populations
  (the existing `rotation_residuals` parity tests continue to pin the
  dispatcher — both arms change together, so they stay green).
- `null9_from_8rows`: for seeded random rank-8 designs, `‖A·v‖ ≈ 0` and unit
  norm, and the span agrees with the eigen path (`|v_elim · v_eig| ≈ 1`);
  the zero design is `None`; a rank-deficient design returns a valid null
  vector. Existing `homography_4pt` exact-reconstruction and degeneracy tests
  stay green.
- Existing focal-vote / estimate-intrinsics suites: values pinned to the old
  arithmetic that legitimately move in last bits get updated pins; every such
  update is listed in the report.

## Bench protocol

The observation dumps from the SIMD pass (session scratchpad,
`simd-profile/{kerry_park,dino_dog_toy}.bin`; format `u64 n, u32 w, u32 h,
u32 ci[n], u32 ii[n], f64 xy[2n]`, little-endian). A scratch example (deleted
before finishing) times `focal_vote_with_options` for `[Pinhole]` and
`[Pinhole, EquidistantFisheye]`, min over 3 warm passes, at
`RAYON_NUM_THREADS=1` and full threads, new kernels vs both diagnostic flags
set. Report the table; the experiment's numbers above are the expectation.

## Spec updates in the same change

- [`focal-vote.md`](../core/geometry/focal-vote.md): the "Vectorized residual
  loops" section becomes the determinism discipline of the kernel —
  lane-per-point vectorization, the shared-table polynomial `acos` with its
  scalar twin, the elimination minimal solvers, platform determinism of those
  paths, and the three environment flags.
- `geometry::simd` module doc: the "no vectorized transcendentals" clause is
  rewritten — transcendentals are polynomial with a scalar twin of identical
  operation order, never libm in one arm and approximation in the other.
