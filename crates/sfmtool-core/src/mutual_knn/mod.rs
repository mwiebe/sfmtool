// Copyright The SfM Tool Authors
// SPDX-License-Identifier: Apache-2.0

//! Mutual-kNN cross-image matcher.
//!
//! A lighter-weight alternative to the background-floor cluster matcher (see
//! [`crate::cluster_match`]). Instead of a per-descriptor adaptive radius and a
//! hard cluster partition, it keeps each descriptor's `k` nearest *cross-image*
//! neighbours and emits the **mutual** ones: an edge `a–b` survives iff `b` is
//! among `a`'s top-`k` and `a` is among `b`'s top-`k`. Mutuality replaces the
//! radius as the precision mechanism — the matcher deliberately keeps more (and
//! wider-baseline) candidate matches than the background floor and leans on the
//! downstream geometric verification (RANSAC two-view geometry) to reject the
//! false ones.
//!
//! An optional **triangle filter** keeps a mutual edge only when at least
//! `triangle_min` third descriptors are mutually matched to *both* endpoints (a
//! 3-clique in the mutual graph) — a cheap, geometry-free multi-view-consistency
//! pre-filter that sheds the dangling/uncorroborated edges that otherwise
//! dominate verification cost.
//!
//! The design, the recall/precision analysis, and the rationale for the
//! defaults live in `specs/core/mutual-knn-matching.md`.

use std::borrow::Cow;
use std::time::Instant;

use ndarray::{Array1, Array2, ArrayView1, ArrayView2};
use rayon::prelude::*;

use crate::cluster_match::{ClusterMatchError, PairMatches};
use crate::kdforest::{KdForestParams, KdForestU8};

#[cfg(test)]
mod tests;

/// Parameters for the mutual-kNN matcher.
#[derive(Clone, Debug)]
pub struct MutualKnnParams {
    /// Nearest cross-image neighbours kept per descriptor. Larger `k` recovers
    /// more wide-baseline matches (higher recall) at the cost of more candidate
    /// matches to verify. Default 12.
    pub k: usize,
    /// Triangle filter: keep a mutual edge only if at least this many third
    /// descriptors are mutually matched to *both* endpoints. `0` disables it.
    /// Default 0.
    pub triangle_min: usize,
    /// kd-forest index build + per-query search budget.
    /// Default [`KdForestParams::accurate`].
    pub forest: KdForestParams,
}

impl Default for MutualKnnParams {
    fn default() -> Self {
        Self {
            k: 12,
            triangle_min: 0,
            forest: KdForestParams::accurate(),
        }
    }
}

/// Euclidean L2 distance between two `u8` descriptor rows (integer squared-L2
/// accumulation, `sqrt` only when reporting — like the forest itself).
fn l2_distance(a: ArrayView1<'_, u8>, b: ArrayView1<'_, u8>) -> f32 {
    let mut acc: i64 = 0;
    for (&x, &y) in a.iter().zip(b.iter()) {
        let diff = x as i64 - y as i64;
        acc += diff * diff;
    }
    (acc as f64).sqrt() as f32
}

/// Number of shared elements of two ascending-sorted slices.
fn sorted_intersection_count(a: &[u32], b: &[u32]) -> usize {
    let (mut i, mut j, mut count) = (0usize, 0usize, 0usize);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                count += 1;
                i += 1;
                j += 1;
            }
        }
    }
    count
}

/// Mutual-kNN cross-image matcher.
///
/// `descriptors` is the `(N, D)` uint8 corpus of every image's SIFT descriptors
/// concatenated image by image; `image_starts` is the CSR offset array (length
/// `n_images + 1`) mapping rows to images, exactly as for
/// [`crate::cluster_match::background_floor_clusters`]. Returns the cross-image
/// matches bucketed by image pair, ready for the `.matches` writer.
pub fn mutual_knn_matches(
    descriptors: ArrayView2<'_, u8>,
    image_starts: &[u32],
    params: &MutualKnnParams,
) -> Result<PairMatches, ClusterMatchError> {
    let n = descriptors.nrows();
    let dim = descriptors.ncols();
    if n == 0 {
        return Err(ClusterMatchError::EmptyCorpus);
    }
    let offsets_valid = image_starts.len() >= 2
        && image_starts[0] == 0
        && image_starts.windows(2).all(|w| w[0] <= w[1])
        && *image_starts.last().unwrap() as usize == n;
    if !offsets_valid {
        return Err(ClusterMatchError::BadOffsets { n });
    }
    let n_images = image_starts.len() - 1;
    let k = params.k.max(1);

    // Row -> owning image. Because the corpus is concatenated image by image,
    // global row order follows image order: for rows a < b in *different*
    // images, image_of[a] < image_of[b].
    let mut image_of = vec![0u32; n];
    for img in 0..n_images {
        let lo = image_starts[img] as usize;
        let hi = image_starts[img + 1] as usize;
        image_of[lo..hi].fill(img as u32);
    }

    let corpus: Cow<'_, [u8]> = match descriptors.as_slice() {
        Some(s) => Cow::Borrowed(s),
        None => Cow::Owned(descriptors.iter().copied().collect()),
    };

    // Phase timing for the perf sweep: set MUTUAL_KNN_PROFILE to log each
    // phase's wall time to stderr. Off (and free) otherwise.
    let prof = std::env::var_os("MUTUAL_KNN_PROFILE").is_some();
    let mark = |label: &str, t: Instant| {
        if prof {
            eprintln!(
                "[mutual-knn] {label:>11}: {:>9.1} ms",
                t.elapsed().as_secs_f64() * 1e3
            );
        }
    };

    // The forest returns nearest neighbours over the whole corpus, including
    // self and same-image hits. Over-query so that, after dropping those, we
    // still have `k` genuine cross-image neighbours for most descriptors
    // (repeated texture inside one image can otherwise crowd them out).
    let query_k = (2 * k + 1).min(n);
    let t = Instant::now();
    let forest = KdForestU8::build(&corpus, n, dim, params.forest);
    mark("build", t);
    let t = Instant::now();
    let (idx, _dist) = forest.search_batch_with_distances(
        &corpus,
        n,
        query_k,
        params.forest.max_leaf_checks,
        None,
    );
    mark("query", t);

    // Per-descriptor top-k cross-image neighbour rows (nearest first), packed
    // into a flat n*k array with u32::MAX padding.
    let t = Instant::now();
    let mut nbr = vec![u32::MAX; n * k];
    nbr.par_chunks_mut(k).enumerate().for_each(|(i, row)| {
        let mut m = 0;
        for c in 0..query_k {
            if m >= k {
                break;
            }
            let j = idx[i * query_k + c];
            if j == u32::MAX || j as usize == i || image_of[j as usize] == image_of[i] {
                continue;
            }
            row[m] = j;
            m += 1;
        }
    });
    mark("neighbours", t);

    // Mutual edges: a–b kept iff b is in a's row and a is in b's row. Dedup with
    // a < b; membership is a short linear scan of the length-k neighbour row.
    let t = Instant::now();
    let mut edges: Vec<(u32, u32)> = (0..n)
        .into_par_iter()
        .flat_map_iter(|a| {
            let arow = &nbr[a * k..a * k + k];
            let mut out: Vec<(u32, u32)> = Vec::new();
            for &b in arow {
                if b == u32::MAX || (a as u32) >= b {
                    continue;
                }
                let brow = &nbr[b as usize * k..b as usize * k + k];
                if brow.contains(&(a as u32)) {
                    out.push((a as u32, b));
                }
            }
            out
        })
        .collect();
    mark("mutual", t);

    // Optional triangle filter: keep a–b only if it closes at least
    // `triangle_min` triangles (a third descriptor mutually matched to both).
    let t = Instant::now();
    if params.triangle_min > 0 {
        let mut adj: Vec<Vec<u32>> = vec![Vec::new(); n];
        for &(a, b) in &edges {
            adj[a as usize].push(b);
            adj[b as usize].push(a);
        }
        adj.par_iter_mut().for_each(|v| v.sort_unstable());
        edges = edges
            .into_par_iter()
            .filter(|&(a, b)| {
                sorted_intersection_count(&adj[a as usize], &adj[b as usize]) >= params.triangle_min
            })
            .collect();
    }
    mark("triangle", t);

    // Expand to image-pair-bucketed matches. For an edge a < b in different
    // images, image_of[a] < image_of[b], so the pair is already ordered.
    let t = Instant::now();
    let mut tuples: Vec<(u32, u32, u32, u32, f32)> = edges
        .par_iter()
        .map(|&(a, b)| {
            let img_lo = image_of[a as usize];
            let img_hi = image_of[b as usize];
            let feat_lo = a - image_starts[img_lo as usize];
            let feat_hi = b - image_starts[img_hi as usize];
            let dist = l2_distance(descriptors.row(a as usize), descriptors.row(b as usize));
            (img_lo, img_hi, feat_lo, feat_hi, dist)
        })
        .collect();
    tuples.par_sort_unstable_by(|x, y| (x.0, x.1, x.2, x.3).cmp(&(y.0, y.1, y.2, y.3)));

    let mut image_index_pairs: Vec<u32> = Vec::new();
    let mut match_counts: Vec<u32> = Vec::new();
    let mut match_feature_indexes: Vec<u32> = Vec::with_capacity(tuples.len() * 2);
    let mut match_descriptor_distances: Vec<f32> = Vec::with_capacity(tuples.len());
    for &(img_lo, img_hi, feat_lo, feat_hi, dist) in &tuples {
        let is_new_pair = match image_index_pairs.rchunks(2).next() {
            Some(last) => last != [img_lo, img_hi],
            None => true,
        };
        if is_new_pair {
            image_index_pairs.extend([img_lo, img_hi]);
            match_counts.push(0);
        }
        *match_counts.last_mut().unwrap() += 1;
        match_feature_indexes.extend([feat_lo, feat_hi]);
        match_descriptor_distances.push(dist);
    }
    mark("expand", t);

    let pair_count = match_counts.len();
    let match_count = match_descriptor_distances.len();
    if prof {
        eprintln!(
            "[mutual-knn]      N={n} k={k} query_k={query_k} pairs={pair_count} matches={match_count}"
        );
    }
    Ok(PairMatches {
        image_index_pairs: Array2::from_shape_vec((pair_count, 2), image_index_pairs).unwrap(),
        match_counts: Array1::from_vec(match_counts),
        match_feature_indexes: Array2::from_shape_vec((match_count, 2), match_feature_indexes)
            .unwrap(),
        match_descriptor_distances: Array1::from_vec(match_descriptor_distances),
    })
}
