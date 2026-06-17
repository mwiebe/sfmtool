# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Mutual-kNN cross-image matching.

Like :mod:`._cluster_matching`, this concatenates every image's SIFT
descriptors into one corpus and matches it directly rather than image pair by
image pair. But instead of a per-descriptor background-floor radius and a hard
cluster partition, it keeps each descriptor's ``k`` nearest *cross-image*
neighbours and emits the **mutual** ones: an edge ``a-b`` survives iff ``b`` is
among ``a``'s top-``k`` and ``a`` is among ``b``'s. Mutuality (not a radius) is
the precision mechanism, so it keeps more — and wider-baseline — candidate
matches and leans on the downstream geometric verification to reject the false
ones. An optional triangle filter keeps only edges corroborated by a third
mutually-matched descriptor.

The kd-forest build, k-NN query, and mutual/triangle logic happen in Rust; see
``specs/core/mutual-knn-matching.md``.
"""

from pathlib import Path
from typing import Optional

import numpy as np

from .._sfmtool import mutual_knn_matches as _rust_mutual_knn_matches
from ..sift.file import SiftReader
from ._cluster_matching import PairArrays


def mutual_knn_match(
    image_paths: list[Path],
    sift_paths: list[Path],
    *,
    k: int = 12,
    triangle_min: int = 0,
    preset: str = "accurate",
    max_feature_count: Optional[int] = None,
) -> PairArrays:
    """Run the mutual-kNN matcher over every image's SIFT descriptors.

    Loads each image's descriptors (capped at ``max_feature_count`` to match the
    feature indices used downstream), concatenates them into one ``(N, 128)``
    uint8 corpus with a CSR ``image_starts`` array, and calls
    ``sfmtool.mutual_knn_matches``. Returns the four parallel pair arrays
    (``image_index_pairs``, ``match_counts``, ``match_feature_indexes``,
    ``match_descriptor_distances``) the ``.matches`` writer consumes — the same
    shape :func:`._cluster_matching.cluster_match` produces.
    """
    assert len(image_paths) == len(sift_paths)

    descriptors = []
    for sift_path in sift_paths:
        with SiftReader(sift_path) as reader:
            descriptors.append(reader.read_descriptors(count=max_feature_count))

    image_starts = np.zeros(len(descriptors) + 1, dtype=np.uint32)
    image_starts[1:] = np.cumsum([len(desc) for desc in descriptors])
    if descriptors:
        corpus = np.ascontiguousarray(np.concatenate(descriptors, axis=0))
    else:
        corpus = np.zeros((0, 128), dtype=np.uint8)

    return PairArrays(
        *_rust_mutual_knn_matches(
            corpus,
            image_starts,
            k=k,
            triangle_min=triangle_min,
            preset=preset,
        )
    )
