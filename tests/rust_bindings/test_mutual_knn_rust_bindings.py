# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Tests for the mutual-kNN cross-image matcher Rust bindings."""

import numpy as np
import pytest

from sfmtool._sfmtool import mutual_knn_matches

N_IMAGES = 4
N_POINTS = 5
N_BACKGROUND = 20
DIM = 128


def _planted_corpus(seed=42):
    """Corpus with N_POINTS planted cross-image points plus background.

    Each image's rows start with the planted observations (base descriptor +
    small jitter), so a planted row's feature index equals its point id;
    N_BACKGROUND random rows per image follow.
    """
    rng = np.random.default_rng(seed)
    bases = rng.integers(0, 256, size=(N_POINTS, DIM), dtype=np.int16)

    blocks = []
    image_starts = [0]
    for _ in range(N_IMAGES):
        jitter = rng.integers(-2, 3, size=(N_POINTS, DIM), dtype=np.int16)
        planted = np.clip(bases + jitter, 0, 255).astype(np.uint8)
        background = rng.integers(0, 256, size=(N_BACKGROUND, DIM), dtype=np.uint8)
        blocks.append(np.vstack([planted, background]))
        image_starts.append(image_starts[-1] + len(blocks[-1]))

    corpus = np.ascontiguousarray(np.vstack(blocks))
    return corpus, np.asarray(image_starts, dtype=np.uint32)


def _exact_kwargs(n):
    """Single-leaf forest config: the k-NN search is exact by construction."""
    return dict(num_trees=1, leaf_size=n, max_leaf_checks=n)


def _candidate_set(pairs, counts, feat_idx):
    """All emitted matches as a set of (img_i, img_j, feat_i, feat_j) tuples."""
    out = set()
    offset = 0
    for k in range(len(pairs)):
        count = int(counts[k])
        for r in range(offset, offset + count):
            out.add(
                (
                    int(pairs[k, 0]),
                    int(pairs[k, 1]),
                    int(feat_idx[r, 0]),
                    int(feat_idx[r, 1]),
                )
            )
        offset += count
    return out


def _planted_matches():
    """Every (image pair, planted point) cross-image match."""
    out = set()
    for i in range(N_IMAGES):
        for j in range(i + 1, N_IMAGES):
            for p in range(N_POINTS):
                out.add((i, j, p, p))
    return out


class TestMutualKnnMatches:
    def test_shapes_dtypes_and_pair_ordering(self):
        corpus, image_starts = _planted_corpus()
        pairs, counts, feat_idx, distances = mutual_knn_matches(
            corpus, image_starts, k=N_IMAGES, **_exact_kwargs(len(corpus))
        )
        assert pairs.dtype == np.uint32 and pairs.shape[1] == 2
        assert counts.dtype == np.uint32
        assert feat_idx.dtype == np.uint32 and feat_idx.shape[1] == 2
        assert distances.dtype == np.float32
        assert counts.sum() == len(feat_idx) == len(distances)
        # Pairs are i < j and sorted ascending by (i, j).
        assert np.all(pairs[:, 0] < pairs[:, 1])
        order = np.lexsort((pairs[:, 1], pairs[:, 0]))
        assert np.array_equal(order, np.arange(len(pairs)))

    def test_recovers_the_planted_cross_image_matches(self):
        corpus, image_starts = _planted_corpus()
        pairs, counts, feat_idx, _ = mutual_knn_matches(
            corpus, image_starts, k=N_IMAGES, **_exact_kwargs(len(corpus))
        )
        got = _candidate_set(pairs, counts, feat_idx)
        # Every planted observation is mutually nearest across images, so each
        # image pair carries its (p, p) match.
        assert _planted_matches() <= got

    def test_distances_are_euclidean_l2(self):
        corpus, image_starts = _planted_corpus()
        pairs, counts, feat_idx, distances = mutual_knn_matches(
            corpus, image_starts, k=N_IMAGES, **_exact_kwargs(len(corpus))
        )
        offset = 0
        for k in range(len(pairs)):
            row_i = int(image_starts[pairs[k, 0]]) + int(feat_idx[offset, 0])
            row_j = int(image_starts[pairs[k, 1]]) + int(feat_idx[offset, 1])
            expected = np.linalg.norm(
                corpus[row_i].astype(np.float64) - corpus[row_j].astype(np.float64)
            )
            assert distances[offset] == pytest.approx(expected, rel=1e-5)
            offset += int(counts[k])

    def test_triangle_filter_keeps_planted_and_does_not_grow_candidates(self):
        corpus, image_starts = _planted_corpus()
        base = mutual_knn_matches(
            corpus,
            image_starts,
            k=N_IMAGES,
            triangle_min=0,
            **_exact_kwargs(len(corpus)),
        )
        filtered = mutual_knn_matches(
            corpus,
            image_starts,
            k=N_IMAGES,
            triangle_min=1,
            **_exact_kwargs(len(corpus)),
        )
        # Planted points span N_IMAGES (>= 3) images, so every planted edge closes
        # triangles and survives the filter.
        assert _planted_matches() <= _candidate_set(
            filtered[0], filtered[1], filtered[2]
        )
        # The filter only ever removes edges.
        assert filtered[1].sum() <= base[1].sum()

    def test_determinism(self):
        corpus, image_starts = _planted_corpus(seed=7)
        a = mutual_knn_matches(corpus, image_starts, k=12, preset="accurate")
        b = mutual_knn_matches(corpus, image_starts, k=12, preset="accurate")
        assert np.array_equal(a[0], b[0])
        assert np.array_equal(a[1], b[1])
        assert np.array_equal(a[2], b[2])

    def test_invalid_inputs_raise(self):
        corpus, image_starts = _planted_corpus()
        with pytest.raises(ValueError):
            mutual_knn_matches(corpus, image_starts, k=0)
        with pytest.raises(ValueError):
            # image_starts that does not end at N.
            mutual_knn_matches(corpus, np.asarray([0, 1], dtype=np.uint32), k=4)
