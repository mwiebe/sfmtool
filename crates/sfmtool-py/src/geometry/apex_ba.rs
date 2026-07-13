// Copyright The SfM Tool Authors
// SPDX-License-Identifier: Apache-2.0

//! EVALUATION binding (branch `apex-ba-eval`): robust bundle adjustment via
//! the `apex-solver` crate, mirroring the scipy `least_squares` pass in
//! `scripts/exp_pinhole_bootstrap.py::solve_round`.
//!
//! Problem: minimize a robust reprojection error over per-image SE3 poses
//! (`x_cam = R·X + t`, +z forward), per-cluster 3D points, and one shared
//! focal `f` (optionally fixed). Observations are centered pixels, so the
//! principal point is identically 0 and never optimized.
//!
//! Expressibility notes (primary evaluation findings):
//!
//! - **Fixed f** is mechanically expressible with the built-in
//!   `ProjectionFactor<PinholeCamera, BundleAdjustment>`: intrinsics are
//!   plain constants (`fx = fy = f`, `cx = cy = 0`), one factor per
//!   observation connecting a pose variable and a point variable. It was
//!   tried first for `opt_f=False`, but its behind-camera convention (zero
//!   residual AND zero Jacobian, Ceres-style) is unusable for our
//!   anchor-free BA: reflecting the entire scene behind the cameras zeroes
//!   almost every residual, so the mirror solution is a near-costless
//!   minimum — on the seoul dataset the final BA collapsed into exactly
//!   that (camera rotation errors ~178 deg vs the reference). Both paths
//!   therefore use the custom factor below, which keeps a large penalty
//!   behind the camera by projecting at a floored depth (`z >= 1e-3 * f`,
//!   the caller's own trim criterion; see `MIN_DEPTH_REL` for why scipy's
//!   exact `max(z, 1e-6)` clamp is numerically fatal to the explicit Schur
//!   backend).
//! - **Free single f** is NOT expressible with the built-ins: the pinhole
//!   intrinsics variable is all-or-nothing 4-parameter `[fx, fy, cx, cy]`
//!   (`INTRINSIC_DIM = 4`), with no way to tie `fx = fy` to one parameter.
//!   `Problem::fix_variable` can pin individual indexes, but it only zeroes
//!   the applied tangent step after the linear solve (the Jacobian columns
//!   stay in the system), i.e. a projected-step approximation rather than
//!   true gating — and it cannot tie parameters at all. `opt_f=True`
//!   therefore uses [`SingleFocalProjection`], a custom `Factor` with the
//!   exact 1-parameter-focal residual and analytic Jacobians (2x6 pose in
//!   the SE3 right-perturbation tangent, 2x3 point, 2x1 focal).
//! - **Loss**: apex has no soft_l1/pseudo-Huber. `BarronGeneralLoss` at
//!   alpha=1 is nominally Charbonnier, but in 1.3.0 its rho/rho'/rho'' are
//!   mutually inconsistent (each implies a different scale), so we implement
//!   scipy's exact soft_l1 as a custom [`LossFunction`] instead — the trait
//!   is public and takes `s = ||r||^2` Ceres-style. One semantic difference
//!   stands: scipy robustifies each residual *component*, apex each 2D
//!   observation block (`s = du^2 + dv^2`); the block form is kept as-is.
//! - **Backend**: explicit sparse Schur complement
//!   (`SparseSchurComplement` + `SchurVariant::Sparse`). At ~130 cameras /
//!   ~10k points the reduced camera system is only ~780 DOF, while plain
//!   sparse Cholesky would factor the full ~31k-DOF system in which the
//!   shared focal column is dense. Landmark elimination requires the point
//!   variables to be literally named `pt_*` (a name-pattern contract of the
//!   solver); the shared focal, named differently, lands in the reduced
//!   camera block where it belongs. Gauge freedom is handled by LM damping,
//!   like scipy's trust region.

use std::borrow::Cow;
use std::collections::HashMap;

use apex_nalgebra::{DMatrix, DVector, Quaternion, UnitQuaternion, Vector3};
use apex_solver::core::loss_functions::LossFunction;
use apex_solver::core::problem::Problem;
use apex_solver::factors::Factor;
use apex_solver::linalg::sparse::SchurVariant;
use apex_solver::manifold::se3::SE3;
use apex_solver::manifold::{LieGroup, ManifoldType};
use apex_solver::{JacobianMode, LevenbergMarquardt, LevenbergMarquardtConfig, LinearSolverType};
use numpy::{PyArray2, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::matching::cluster::extract_u32_1d;

/// Relative depth floor: `z` is clamped to `MIN_DEPTH_REL * f` (with a tiny
/// absolute backstop) and the smooth projection formulas are evaluated at
/// the floored depth. `1e-3 * f` is the same behind-camera criterion the
/// caller's trim step uses (`depth > 1e-3 * f`).
///
/// Why not scipy's exact `max(z, 1e-6)` with a vanishing `d/dz`? Two
/// numerically fatal interactions with this solver were observed on seoul:
/// residuals of ~1e9 px put ~1e17 entries into the normal equations, and
/// the zero z-column makes the per-landmark 3x3 Hessian blocks rank-2, so
/// the explicit Schur elimination divides by the LM damping and injects
/// ~1e15 (or NaN) entries into the reduced system ("Schur complement
/// singular after 5 regularization attempts"). The floored-depth form keeps
/// a strong penalty against behind-camera configurations (hundreds-of-px
/// residuals, exactly what prevents the whole-scene reflection collapse)
/// while keeping every block full-rank and bounded.
const MIN_DEPTH_REL: f64 = 1e-3;

/// scipy's `soft_l1` loss, Ceres-style on the squared block norm `s`:
/// `rho(s) = 2 c^2 (sqrt(1 + s/c^2) - 1)` with `rho'(0) = 1`.
struct SoftL1Loss {
    c2: f64,
}

impl LossFunction for SoftL1Loss {
    #[inline]
    fn evaluate(&self, s: f64) -> [f64; 3] {
        let t = 1.0 + s / self.c2;
        let sq = t.sqrt();
        [
            2.0 * self.c2 * (sq - 1.0),
            1.0 / sq,
            -1.0 / (2.0 * self.c2 * sq * sq * sq),
        ]
    }
}

/// Custom projection factor with a single shared focal parameter.
///
/// Residual: `(f·x/z - u, f·y/z - v)` for `p_cam = (x, y, z) = R·X + t`.
/// Connected variables: `[pose (SE3), point (R3)]` plus `[focal (R1)]` when
/// the focal is free. Jacobians are analytic; the pose block uses the
/// right-perturbation convention of `apex_manifolds::se3` (`d p_cam/d rho =
/// R`, `d p_cam/d theta = -R·[X]x`), identical to the built-in factor.
struct SingleFocalProjection {
    obs: [f64; 2],
    /// `Some(f)` bakes the focal in as a constant (no focal variable).
    fixed_f: Option<f64>,
}

impl Factor for SingleFocalProjection {
    fn linearize(
        &self,
        params: &[DVector<f64>],
        compute_jacobian: bool,
    ) -> (DVector<f64>, Option<DMatrix<f64>>) {
        let pose = SE3::from(params[0].clone());
        let pt = Vector3::new(params[1][0], params[1][1], params[1][2]);
        let f = self.fixed_f.unwrap_or_else(|| params[2][0]);
        let ncols = if self.fixed_f.is_some() { 9 } else { 10 };

        let p_cam = pose.act(&pt, None, None);
        // Behind-camera observations keep a large (hundreds of px) residual
        // rather than the Ceres zero-residual convention — with zero
        // residuals, reflecting the whole gauge-free scene behind the
        // cameras is a near-costless minimum and the built-in
        // ProjectionFactor was observed to collapse into exactly that. The
        // smooth formulas are evaluated at the floored depth (full-rank
        // Jacobians; see MIN_DEPTH_REL for why d/dz must not vanish).
        let inv_z = 1.0 / p_cam.z.max((MIN_DEPTH_REL * f).max(1e-6));
        let xn = p_cam.x * inv_z;
        let yn = p_cam.y * inv_z;
        let residual = DVector::from_vec(vec![f * xn - self.obs[0], f * yn - self.obs[1]]);

        let jacobian = compute_jacobian.then(|| {
            // d(u,v)/d(p_cam), 2x3, at the floored depth.
            let g = f * inv_z;
            let d_uv_d_pcam =
                apex_nalgebra::SMatrix::<f64, 2, 3>::new(g, 0.0, -g * xn, 0.0, g, -g * yn);
            let rot = pose.rotation_so3().rotation_matrix();
            // Point block; also the translation half of the pose block.
            let d_uv_d_point = d_uv_d_pcam * rot;
            // Rotation half: d p_cam/d theta = -R·[X]x.
            let skew = apex_nalgebra::Matrix3::new(
                0.0, -pt.z, pt.y, //
                pt.z, 0.0, -pt.x, //
                -pt.y, pt.x, 0.0,
            );
            let d_uv_d_theta = d_uv_d_pcam * (-rot * skew);

            let mut jac = DMatrix::zeros(2, ncols);
            for r in 0..2 {
                for c in 0..3 {
                    jac[(r, c)] = d_uv_d_point[(r, c)]; // pose: rho
                    jac[(r, 3 + c)] = d_uv_d_theta[(r, c)]; // pose: theta
                    jac[(r, 6 + c)] = d_uv_d_point[(r, c)]; // point
                }
            }
            if self.fixed_f.is_none() {
                jac[(0, 9)] = xn;
                jac[(1, 9)] = yn;
            }
            jac
        });
        (residual, jacobian)
    }

    fn get_dimension(&self) -> usize {
        2
    }
}

/// SE3 initial value `[tx, ty, tz, qw, qx, qy, qz]` from an axis-angle
/// rotation vector and a translation (scipy `Rotation.from_rotvec`
/// convention: both describe the same `x_cam = R·X + t`).
fn se3_value(rv: &[f64], tv: &[f64]) -> DVector<f64> {
    let q = UnitQuaternion::from_scaled_axis(Vector3::new(rv[0], rv[1], rv[2]));
    DVector::from_vec(vec![tv[0], tv[1], tv[2], q.w, q.i, q.j, q.k])
}

/// Back-conversion: `[tx, ty, tz, qw, qx, qy, qz]` to (rotvec, translation),
/// canonicalized to `qw >= 0` so the rotvec magnitude is <= pi (matching
/// scipy's `as_rotvec`).
fn se3_to_rvec_tvec(v: &DVector<f64>) -> ([f64; 3], [f64; 3]) {
    let sign = if v[3] < 0.0 { -1.0 } else { 1.0 };
    let q = UnitQuaternion::from_quaternion(Quaternion::new(
        sign * v[3],
        sign * v[4],
        sign * v[5],
        sign * v[6],
    ));
    let rv = q.scaled_axis();
    ([rv.x, rv.y, rv.z], [v[0], v[1], v[2]])
}

fn check_rows3(name: &str, arr: &PyReadonlyArray2<'_, f64>) -> PyResult<()> {
    if arr.shape()[1] != 3 {
        return Err(PyValueError::new_err(format!(
            "{name} must have shape (N, 3), got (N, {})",
            arr.shape()[1]
        )));
    }
    Ok(())
}

/// One robust bundle-adjustment pass via apex-solver (EVALUATION binding).
///
/// Mirrors the scipy pass in ``exp_pinhole_bootstrap.solve_round``: robust
/// (soft_l1, per-observation) reprojection error over per-image SE3 poses,
/// per-point 3D positions, and one shared focal, solved with
/// Levenberg-Marquardt over an explicit sparse Schur complement.
///
/// Args:
///     obs_points: (K,) uint32 compact point index per observation.
///     obs_images: (K,) uint32 compact image index per observation.
///     obs_xy: (K, 2) float64 centered pixel observations.
///     rvec: (N, 3) float64 axis-angle rotations, ``cam_from_world``.
///     tvec: (N, 3) float64 translations, ``cam_from_world``.
///     points: (P, 3) float64 3D points.
///     f: Shared focal length in pixels.
///     opt_f: Optimize the focal (True) or hold it fixed (False).
///     loss_scale: soft_l1 scale in pixels (scipy ``f_scale``).
///     max_iterations: LM iteration cap (compare scipy ``max_nfev``).
///
/// Returns:
///     Dict with ``f`` (float), ``rvec`` (N, 3), ``tvec`` (N, 3),
///     ``points`` (P, 3), ``iterations`` (int), ``final_cost`` (float).
#[pyfunction]
#[pyo3(signature = (obs_points, obs_images, obs_xy, rvec, tvec, points, f, opt_f, loss_scale, max_iterations))]
#[allow(clippy::too_many_arguments)]
pub fn bundle_adjust_apex(
    py: Python<'_>,
    obs_points: &Bound<'_, PyAny>,
    obs_images: &Bound<'_, PyAny>,
    obs_xy: PyReadonlyArray2<'_, f64>,
    rvec: PyReadonlyArray2<'_, f64>,
    tvec: PyReadonlyArray2<'_, f64>,
    points: PyReadonlyArray2<'_, f64>,
    f: f64,
    opt_f: bool,
    loss_scale: f64,
    max_iterations: usize,
) -> PyResult<Py<PyAny>> {
    let obs_points = extract_u32_1d(obs_points, "obs_points")?;
    let obs_images = extract_u32_1d(obs_images, "obs_images")?;
    let op: Cow<'_, [u32]> = to_contiguous!(obs_points);
    let oi: Cow<'_, [u32]> = to_contiguous!(obs_images);
    if obs_xy.shape()[1] != 2 {
        return Err(PyValueError::new_err(format!(
            "obs_xy must have shape (K, 2), got (K, {})",
            obs_xy.shape()[1]
        )));
    }
    check_rows3("rvec", &rvec)?;
    check_rows3("tvec", &tvec)?;
    check_rows3("points", &points)?;
    let k = op.len();
    if oi.len() != k || obs_xy.shape()[0] != k {
        return Err(PyValueError::new_err(
            "obs_points, obs_images and obs_xy must have the same length",
        ));
    }
    let n_img = rvec.shape()[0];
    let n_pt = points.shape()[0];
    if tvec.shape()[0] != n_img {
        return Err(PyValueError::new_err("rvec and tvec must have equal rows"));
    }
    if !loss_scale.is_finite() || loss_scale <= 0.0 {
        return Err(PyValueError::new_err("loss_scale must be positive"));
    }
    if op.iter().any(|&p| p as usize >= n_pt) || oi.iter().any(|&i| i as usize >= n_img) {
        return Err(PyValueError::new_err(
            "observation indexes exceed the pose/point array sizes",
        ));
    }
    let xy: Cow<'_, [f64]> = to_contiguous!(obs_xy);
    let rv: Cow<'_, [f64]> = to_contiguous!(rvec);
    let tv: Cow<'_, [f64]> = to_contiguous!(tvec);
    let pt: Cow<'_, [f64]> = to_contiguous!(points);

    // Variable names. Sorted-name order defines the Jacobian column layout,
    // and the `pt_` prefix is the solver's Schur landmark-elimination
    // contract; zero-padding keeps everything deterministic.
    let cam_names: Vec<String> = (0..n_img).map(|i| format!("cam_{i:06}")).collect();
    let pt_names: Vec<String> = (0..n_pt).map(|p| format!("pt_{p:06}")).collect();
    const FOCAL: &str = "focal";

    let mut problem = Problem::new(JacobianMode::Sparse);
    let c2 = loss_scale * loss_scale;
    for j in 0..k {
        let cam = cam_names[oi[j] as usize].as_str();
        let ptn = pt_names[op[j] as usize].as_str();
        let obs = [xy[2 * j], xy[2 * j + 1]];
        let loss: Option<Box<dyn LossFunction + Send>> = Some(Box::new(SoftL1Loss { c2 }));
        // The custom factor is used for BOTH the free- and fixed-focal
        // paths: the built-in ProjectionFactor is mechanically expressible
        // for a fixed focal but its zero-residual cheirality convention
        // collapses gauge-free BAs (see the module docs).
        if opt_f {
            let factor = SingleFocalProjection { obs, fixed_f: None };
            problem.add_residual_block(&[cam, ptn, FOCAL], Box::new(factor), loss);
        } else {
            let factor = SingleFocalProjection {
                obs,
                fixed_f: Some(f),
            };
            problem.add_residual_block(&[cam, ptn], Box::new(factor), loss);
        }
    }

    let mut initial: HashMap<String, (ManifoldType, DVector<f64>)> = HashMap::new();
    for i in 0..n_img {
        initial.insert(
            cam_names[i].clone(),
            (
                ManifoldType::SE3,
                se3_value(&rv[3 * i..3 * i + 3], &tv[3 * i..3 * i + 3]),
            ),
        );
    }
    for p in 0..n_pt {
        initial.insert(
            pt_names[p].clone(),
            (
                ManifoldType::RN,
                DVector::from_vec(pt[3 * p..3 * p + 3].to_vec()),
            ),
        );
    }
    if opt_f {
        initial.insert(
            FOCAL.to_string(),
            (ManifoldType::RN, DVector::from_vec(vec![f])),
        );
    }

    let config = LevenbergMarquardtConfig::new()
        .with_linear_solver_type(LinearSolverType::SparseSchurComplement)
        .with_schur_variant(SchurVariant::Sparse)
        .with_max_iterations(max_iterations)
        .with_cost_tolerance(1e-8)
        .with_parameter_tolerance(1e-8);

    let result = py
        .detach(|| {
            let mut solver = LevenbergMarquardt::with_config(config);
            solver.optimize(&problem, &initial)
        })
        .map_err(|e| PyRuntimeError::new_err(format!("apex-solver failed: {e}")))?;

    let vars = &result.parameters;
    let missing = |name: &str| PyRuntimeError::new_err(format!("variable {name} missing"));
    let mut rvec_out: Vec<Vec<f64>> = Vec::with_capacity(n_img);
    let mut tvec_out: Vec<Vec<f64>> = Vec::with_capacity(n_img);
    for name in &cam_names {
        let v = vars.get(name).ok_or_else(|| missing(name))?.to_vector();
        let (r, t) = se3_to_rvec_tvec(&v);
        rvec_out.push(r.to_vec());
        tvec_out.push(t.to_vec());
    }
    let mut pts_out: Vec<Vec<f64>> = Vec::with_capacity(n_pt);
    for name in &pt_names {
        let v = vars.get(name).ok_or_else(|| missing(name))?.to_vector();
        pts_out.push(vec![v[0], v[1], v[2]]);
    }
    let f_out = if opt_f {
        vars.get(FOCAL).ok_or_else(|| missing(FOCAL))?.to_vector()[0]
    } else {
        f
    };

    let err = |e: numpy::FromVecError| PyValueError::new_err(e.to_string());
    let dict = PyDict::new(py);
    dict.set_item("f", f_out)?;
    dict.set_item("rvec", PyArray2::from_vec2(py, &rvec_out).map_err(err)?)?;
    dict.set_item("tvec", PyArray2::from_vec2(py, &tvec_out).map_err(err)?)?;
    dict.set_item("points", PyArray2::from_vec2(py, &pts_out).map_err(err)?)?;
    dict.set_item("iterations", result.iterations)?;
    dict.set_item("final_cost", result.final_cost)?;
    Ok(dict.into_any().unbind())
}

pub fn register(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(pyo3::wrap_pyfunction!(bundle_adjust_apex, m)?)?;
    Ok(())
}
