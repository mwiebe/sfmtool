# Staged robust bundle adjustment (shared camera)

**Status:** Implemented —
`crates/sfmtool-core/src/geometry/bundle_adjust.rs` (`bundle_adjust`,
`BaSchedule`, `BundleAdjustment`; tests in `bundle_adjust/tests.rs`), PyO3
binding in `crates/sfmtool-py/src/geometry/bundle_adjust.rs`
(`sfmtool._sfmtool.geometry.bundle_adjust`), Python tests in
`tests/rust_bindings/test_bundle_adjust_rust_bindings.py`.

## Purpose

The staged robust bundle adjustment used by the cluster pinhole bootstrap
(`specs/core/cluster-pinhole-bootstrap.md`,
`scripts/exp_fast_pinhole.py` / `scripts/exp_pinhole_bootstrap.py`): given
images sharing one camera model, camera poses, world points, and pixel
observations tying them together, jointly refine the poses and points (and
optionally the shared focal length) by minimizing robust pixel reprojection
error over a trim schedule with inter-round retriangulation.

This is the optimizer that the trimmed pose-only refinement
(`crates/sfmtool-core/src/geometry/pose_refine.rs`) is the single-pose
special case of. It replaces the experiment scripts'
`scipy.optimize.least_squares` BA, whose Python-side residual and sparsity
handling dominated the bootstrap's wall-clock.

## Definitions

- `n_img` **images** sharing one `CameraIntrinsics`, each with a
  world-to-camera pose `(R_i, t_i)` in the canonical convention
  (`x_cam = R·X + t`; the camera looks along `−Z`, a point in front has
  `z < 0`), rotations supplied as WXYZ unit quaternions.
- `n_pt` world **points** `X_p` (canonical world frame). Points may be
  non-finite (`NaN`) — their observations are invalid until a
  retriangulation round replaces them.
- `n_obs` **observations** `(image, point, uv)` with `uv` the observed full
  (un-centered) pixel position.
- A **track** is the set of observations of one point.

The state arrays are full-sized: images and points never touched by an
observation pass through unchanged (the solve compacts internally over what
the observations reference).

## The staged loop

```rust
pub struct BaSchedule {
    pub trim_px: f64,     // pre-round trim threshold on the residual norm
    pub loss_scale: f64,  // soft-L1 scale for the round's solve, px
}

pub fn bundle_adjust(
    cam: &CameraIntrinsics,          // shared model; carries the initial focal
    quats: &mut [UnitQuaternion<f64>],   // n_img, world-to-camera
    trans: &mut [Vector3<f64>],          // n_img
    points: &mut [[f64; 3]],             // n_pt (NaN allowed)
    uv: &[[f64; 2]],                     // n_obs
    obs_img: &[u32],                     // n_obs
    obs_pt: &[u32],                      // n_obs
    opt_f: bool,
    schedule: &[BaSchedule],             // default 50/5 → 12/2 → 4/1
    max_iters: usize,                    // LM iterations per round
    min_track: usize,                    // trim survivors per point (2)
    min_obs: usize,                      // degenerate-exit floor (12)
) -> BundleAdjustment;                   // { focal, residual_norms }
```

Per schedule round, mirroring the experiment scripts exactly:

1. **Retriangulate (rounds after the first).** Rebuild *every* point from
   *all* supplied observations at the current poses: world rays
   `R_iᵀ · pixel_to_ray(uv)` and centers `−R_iᵀ t_i` per observation,
   grouped by point, through
   [`reconstruction::triangulation::triangulate_batch`]. A track with fewer
   than 2 observations becomes `NaN`; a point with no observations at all
   becomes `NaN` too (the callers refill from their full observation set —
   the "refill after BA" rule of the bootstrap spec). Re-admission is the
   point: observations a bad init lost re-enter once the refined cameras
   explain them.
2. **Trim.** Keep observations with residual norm `< trim_px`, in-front
   depth `> 1e-3 · f` (canonical depth is `−z_cam`), and a finite point;
   then drop observations of points with fewer than `min_track` survivors.
   If fewer than `min_obs` observations survive, return degenerate: state
   passes through, `residual_norms` all `+∞` (the fast bootstrap's
   "wildly wrong focal" guard).
3. **Solve.** One robust sparse Levenberg–Marquardt solve (below) over the
   kept observations at the round's `loss_scale`.

After the last round, `residual_norms` is the unweighted reprojection
residual norm of **every supplied observation** at the final state (`+∞`
where invalid), so callers tally inlier fractions against denominators of
their own choosing.

## The solve

Levenberg–Marquardt over a local parameterization, minimizing the soft-L1
robust cost

```
cost = Σ_k s² · ρ(‖r_k‖² / s²),   ρ(z) = 2·(√(1 + z) − 1),   s = loss_scale
```

- **Parameters.** Per touched image a local `SO(3) × ℝ³` perturbation
  (`R ← exp(δθ)·R`, `t ← t + δt`); per touched point `X ← X + δX`; when
  `opt_f`, the shared focal `f ← f + δf`. Focal optimization requires
  `SIMPLE_PINHOLE` (single focal, no distortion), where
  `∂(u, v)/∂f = ((u − cx)/f, (v − cy)/f)` exactly; other models are
  rejected at the binding.
- **Jacobian.** Analytic throughout: the projection block from
  `CameraIntrinsics::ray_to_pixel_with_jacobian` composed with `−[R·X]ₓ`
  (rotation), `I₃` (translation), and `R` (point) blocks, exactly as in
  `pose_refine.rs`. An observation whose point is behind the camera /
  outside the model domain contributes residual `(1e6, 0)` with a zero
  Jacobian row — penalized, never steering.
- **Robust weighting.** First-order IRLS: residual and Jacobian rows scale
  by `w = ρ'(z)^½ = (1 + z)^(−¼)`; the true robust cost (not the weighted
  surrogate) decides step acceptance.
- **Schur complement.** Points are eliminated: per-point 3×3 blocks are
  inverted directly and the reduced camera system
  (`[f? | 6·n_im]`, dense) is solved by LU; point updates back-substitute.
  Rejected steps re-damp and re-solve from the same linearization (no
  re-evaluation), with Marquardt scaling `λ·diag(JᵀJ)` for the
  `x_scale="jac"` parameter-scale invariance of the scipy original.
- **Termination.** `max_iters` accepted-step budget per round; stop early
  when an accepted step improves the cost by less than `1e-8` relative, or
  when no damping in a bounded ladder (12 doublings) finds a downhill step.

## Bindings

```python
bundle_adjust(
    camera,                    # CameraIntrinsics shared by all images (initial f)
    quaternions_wxyz,          # (n_img, 4) world-to-camera (WXYZ)
    translations,              # (n_img, 3)
    points,                    # (n_pt, 3), NaN allowed
    uv,                        # (n_obs, 2)
    obs_image,                 # (n_obs,) uint32
    obs_point,                 # (n_obs,) uint32
    opt_f=False,               # requires SIMPLE_PINHOLE
    schedule=[(50.0, 5.0), (12.0, 2.0), (4.0, 1.0)],
    max_iters=60,
    min_track=2,
    min_obs=12,
) -> dict                      # focal, quaternions_wxyz (n_img, 4),
                               # translations (n_img, 3), points (n_pt, 3),
                               # residual_norms (n_obs,)
```

Shapes are validated like `reprojection_residuals`; observation indices out
of range raise. The returned arrays are new (inputs are not mutated from
Python's point of view).

## Testing requirements

- **Perfect-data fixpoint**: synthetic poses/points/observations with zero
  noise stay put (cost already ~0, parameters unchanged to tolerance).
- **Noise recovery**: perturbed poses and points recover the ground truth
  to sub-pixel reprojection on synthetic data; with `opt_f`, a focal
  started 20% off converges to the true value.
- **Robustness**: a contaminated fraction of junk observations does not
  pull the solution (soft-L1 + trim schedule), and the junk ends with
  large `residual_norms` while inliers end small.
- **Trim/track semantics**: an observation set where trimming leaves a
  point with one survivor drops that point's observations from the solve;
  fewer than `min_obs` survivors returns the degenerate all-∞ result with
  the state passed through.
- **Retriangulation re-admission**: a `NaN` point with ≥ 2 observations is
  reborn in round 2 and its observations participate thereafter.
- **Pass-through**: images/points not referenced by any observation are
  returned bit-identical.
- **Binding parity**: the Python binding reproduces the Rust result on a
  small synthetic problem (`tests/rust_bindings/`).

## Non-goals

- Per-image or per-observation camera models — one shared
  `CameraIntrinsics`.
- Optimizing distortion or principal point; `opt_f` covers the single
  shared focal only.
- Gauge fixing, covariance estimation, or constraint handling — callers
  own the gauge (the bootstrap's evaluation aligns by similarity anyway).
- Replacing the production solvers (`sfm solve` wraps COLMAP/GLOMAP); this
  kernel serves the bootstrap experiments and whatever grows out of them.
