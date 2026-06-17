// Copyright The SfM Tool Authors
// SPDX-License-Identifier: Apache-2.0

//! Python binding for the mutual-kNN cross-image matcher.
//!
//! Exposes [`sfmtool_core::mutual_knn`] as `mutual_knn_matches`. It returns the
//! same per-image-pair match arrays `clusters_to_pair_matches` does, so the
//! `.matches` pipeline can consume either matcher's output identically. See
//! `specs/core/mutual-knn-matching.md`.

use std::borrow::Cow;

use ndarray::ArrayView2;
use numpy::{PyArrayMethods, PyReadonlyArray1, PyUntypedArrayMethods};
use pyo3::prelude::*;

use sfmtool_core::mutual_knn::{self, MutualKnnParams};

use crate::py_kdforest::{extract_u8_2d, resolve_forest_params};

/// Extract a 1-D `uint32` array, with a clear error if the dtype is wrong.
fn extract_u32_1d<'py>(
    arr: &Bound<'py, PyAny>,
    what: &str,
) -> PyResult<PyReadonlyArray1<'py, u32>> {
    arr.extract::<PyReadonlyArray1<u32>>().map_err(|_| {
        let dtype = arr
            .getattr("dtype")
            .and_then(|d| d.getattr("name"))
            .and_then(|n| n.extract::<String>())
            .unwrap_or_else(|_| "?".to_string());
        pyo3::exceptions::PyTypeError::new_err(format!(
            "{what} must be a 1-D uint32 array, got {dtype}"
        ))
    })
}

/// Mutual-kNN cross-image matcher.
///
/// Args:
///     descriptors: (N, 128) uint8 corpus, every image's SIFT descriptors
///         concatenated image by image.
///     image_starts: (n_images + 1,) uint32 CSR offsets; image i owns rows
///         image_starts[i]:image_starts[i+1].
///     k: Nearest cross-image neighbours kept per descriptor (default 12).
///         Larger k recovers more wide-baseline matches at more candidates to
///         verify.
///     triangle_min: Keep a mutual edge only if at least this many third
///         descriptors are mutually matched to both endpoints (a 3-clique in
///         the mutual graph). 0 disables the filter (default 0).
///     preset / num_trees / leaf_size / max_leaf_checks / seed: forest config,
///         same meaning as KdForest. The default preset is "accurate".
///
/// Returns:
///     Tuple (image_index_pairs, match_counts, match_feature_indexes,
///     match_descriptor_distances) — identical in shape to
///     `clusters_to_pair_matches`:
///     - image_index_pairs: (P, 2) uint32 sorted pairs with i < j.
///     - match_counts: (P,) uint32 matches per pair.
///     - match_feature_indexes: (M, 2) uint32 feature pairs grouped by pair.
///     - match_descriptor_distances: (M,) float32 Euclidean L2 distances.
#[pyfunction]
#[pyo3(signature = (descriptors, image_starts, k=12, triangle_min=0,
                    preset=None, num_trees=None, leaf_size=None,
                    max_leaf_checks=None, seed=None))]
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn mutual_knn_matches(
    py: Python<'_>,
    descriptors: &Bound<'_, PyAny>,
    image_starts: &Bound<'_, PyAny>,
    k: usize,
    triangle_min: usize,
    preset: Option<&str>,
    num_trees: Option<usize>,
    leaf_size: Option<usize>,
    max_leaf_checks: Option<usize>,
    seed: Option<u64>,
) -> PyResult<(Py<PyAny>, Py<PyAny>, Py<PyAny>, Py<PyAny>)> {
    if k == 0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "k must be at least 1",
        ));
    }
    let descriptors = extract_u8_2d(descriptors, "descriptors")?;
    let shape = descriptors.shape();
    let (n, dim) = (shape[0], shape[1]);
    if dim != 128 {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "descriptors must be (N, 128); got width {dim}"
        )));
    }
    let data: Cow<'_, [u8]> = to_contiguous!(descriptors);
    let image_starts = extract_u32_1d(image_starts, "image_starts")?;
    let starts: Cow<'_, [u32]> = to_contiguous!(image_starts);

    let forest = resolve_forest_params(
        preset,
        "accurate",
        num_trees,
        leaf_size,
        max_leaf_checks,
        seed,
    )?;
    let params = MutualKnnParams {
        k,
        triangle_min,
        forest,
    };

    let pairs = py
        .detach(|| {
            let view = ArrayView2::from_shape((n, dim), data.as_ref()).expect("contiguous corpus");
            mutual_knn::mutual_knn_matches(view, &starts, &params)
        })
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    let pair_count = pairs.image_index_pairs.nrows();
    let match_count = pairs.match_feature_indexes.nrows();
    let image_index_pairs =
        numpy::PyArray1::from_vec(py, pairs.image_index_pairs.into_raw_vec_and_offset().0)
            .reshape([pair_count, 2])?;
    let match_counts =
        numpy::PyArray1::from_vec(py, pairs.match_counts.into_raw_vec_and_offset().0);
    let match_feature_indexes =
        numpy::PyArray1::from_vec(py, pairs.match_feature_indexes.into_raw_vec_and_offset().0)
            .reshape([match_count, 2])?;
    let match_descriptor_distances = numpy::PyArray1::from_vec(
        py,
        pairs.match_descriptor_distances.into_raw_vec_and_offset().0,
    );
    Ok((
        image_index_pairs.into_any().unbind(),
        match_counts.into_any().unbind(),
        match_feature_indexes.into_any().unbind(),
        match_descriptor_distances.into_any().unbind(),
    ))
}
