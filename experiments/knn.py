# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Dependency-free exact k-NN via chunked BLAS GEMM.

Python 3.14 (pinned here) is too new for conda-forge faiss / numba-based ANN,
so the experiments compute *exact* Euclidean nearest neighbours with the
standard ``|a-b|^2 = |a|^2 + |b|^2 - 2 a.b`` trick: the ``a.b`` term is one big
matrix multiply that BLAS handles efficiently.  Done in row-chunks so the
distance block never has to be materialised all at once — this scales to the
~340k-descriptor datasets where a literal 128-D KD-tree (and even sklearn's
brute path) bog down.
"""

from __future__ import annotations

import numpy as np


def knn_indices(
    corpus: np.ndarray,
    queries: np.ndarray,
    k: int,
    *,
    batch: int = 1024,
    drop_first: bool = False,
) -> np.ndarray:
    """Indices of the ``k`` nearest corpus rows for each query row.

    Args:
        corpus: ``(N, D)`` float32 array searched against.
        queries: ``(Q, D)`` float32 array of query rows.
        k: number of neighbours to return per query.
        batch: query rows processed per chunk (caps peak memory).
        drop_first: if True, drop each query's closest hit (its own row when
            ``queries`` is a slice of ``corpus``) and return the next ``k``.

    Returns:
        ``(Q, k)`` int64 array of corpus indices, sorted nearest-first.
    """
    corpus = np.ascontiguousarray(corpus, dtype=np.float32)
    queries = np.ascontiguousarray(queries, dtype=np.float32)
    n = len(corpus)
    fetch = min(k + (1 if drop_first else 0), n)

    corpus_sq = np.einsum("ij,ij->i", corpus, corpus)  # |b|^2, shape (N,)
    out = np.empty((len(queries), k), dtype=np.int64)

    for start in range(0, len(queries), batch):
        qb = queries[start : start + batch]
        # squared distances (Q_b, N); |a|^2 is constant per row so it can't
        # change the ordering and is skipped.
        cross = qb @ corpus.T
        d2 = corpus_sq[None, :] - 2.0 * cross
        # partial sort: cheapest `fetch` per row, then order just those.
        part = np.argpartition(d2, fetch - 1, axis=1)[:, :fetch]
        rows = np.arange(len(qb))[:, None]
        order = np.argsort(d2[rows, part], axis=1)
        nearest = part[rows, order]
        out[start : start + batch] = (
            nearest[:, 1 : 1 + k] if drop_first else nearest[:, :k]
        )

    return out
