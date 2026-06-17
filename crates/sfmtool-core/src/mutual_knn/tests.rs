// Copyright The SfM Tool Authors
// SPDX-License-Identifier: Apache-2.0

use super::*;
use crate::cluster_match::ClusterMatchError;
use crate::kdforest::KdForestParams;
use ndarray::Array2;

const D: usize = 128;

fn params(k: usize, triangle_min: usize) -> MutualKnnParams {
    MutualKnnParams {
        k,
        triangle_min,
        forest: KdForestParams::accurate(),
    }
}

fn corpus(rows: &[Vec<u8>]) -> Array2<u8> {
    let n = rows.len();
    let flat: Vec<u8> = rows.iter().flatten().copied().collect();
    Array2::from_shape_vec((n, D), flat).unwrap()
}

/// A track descriptor: the same vector in every image (distance 0 across views).
fn track() -> Vec<u8> {
    vec![0u8; D]
}

/// A "filler" descriptor: a distinct spike far from the track and from the other
/// fillers, so it forms no mutual edge with anything.
fn filler(dim: usize) -> Vec<u8> {
    let mut v = vec![0u8; D];
    v[dim] = 200;
    v
}

#[test]
fn mutual_edges_are_the_track_across_images() {
    // 3 images, 2 features each: feature 0 is the shared track, feature 1 is a
    // per-image filler. Only the track should produce mutual cross-image edges.
    let rows = vec![track(), filler(0), track(), filler(1), track(), filler(2)];
    let c = corpus(&rows);
    let starts = [0u32, 2, 4, 6];

    let pm = mutual_knn_matches(c.view(), &starts, &params(2, 0)).unwrap();

    // Exactly the three cross-image track pairs, ascending by (i, j).
    assert_eq!(
        pm.image_index_pairs,
        Array2::from_shape_vec((3, 2), vec![0, 1, 0, 2, 1, 2]).unwrap()
    );
    assert_eq!(pm.match_counts.to_vec(), vec![1, 1, 1]);
    // Every match is feature 0 (the track) in both images.
    for row in pm.match_feature_indexes.rows() {
        assert_eq!(row[0], 0);
        assert_eq!(row[1], 0);
    }
    // Track features are identical, so the L2 distance is 0.
    assert!(pm.match_descriptor_distances.iter().all(|&d| d == 0.0));
}

#[test]
fn triangle_filter_keeps_a_three_image_clique() {
    // One feature per image, all identical -> a 3-clique in the mutual graph.
    let rows = vec![track(), track(), track()];
    let c = corpus(&rows);
    let starts = [0u32, 1, 2, 3];

    let kept = mutual_knn_matches(c.view(), &starts, &params(2, 1)).unwrap();
    // All three edges close a triangle, so all survive.
    assert_eq!(kept.match_counts.to_vec(), vec![1, 1, 1]);
}

#[test]
fn triangle_filter_drops_a_dangling_two_image_edge() {
    // One feature per image, two images: a mutual edge with no third corner.
    let rows = vec![track(), track()];
    let c = corpus(&rows);
    let starts = [0u32, 1, 2];

    let no_filter = mutual_knn_matches(c.view(), &starts, &params(2, 0)).unwrap();
    assert_eq!(no_filter.match_counts.sum(), 1);

    let filtered = mutual_knn_matches(c.view(), &starts, &params(2, 1)).unwrap();
    assert_eq!(filtered.match_counts.sum(), 0);
    assert_eq!(filtered.image_index_pairs.nrows(), 0);
}

#[test]
fn empty_corpus_is_an_error() {
    let c = Array2::<u8>::zeros((0, D));
    let starts = [0u32];
    assert!(matches!(
        mutual_knn_matches(c.view(), &starts, &params(4, 0)),
        Err(ClusterMatchError::EmptyCorpus)
    ));
}

#[test]
fn bad_offsets_are_an_error() {
    let rows = vec![track(), track()];
    let c = corpus(&rows);
    // Offsets that do not end at N = 2.
    let starts = [0u32, 1];
    assert!(matches!(
        mutual_knn_matches(c.view(), &starts, &params(2, 0)),
        Err(ClusterMatchError::BadOffsets { .. })
    ));
}
