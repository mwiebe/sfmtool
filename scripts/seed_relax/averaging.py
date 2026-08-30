# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Camera centres from pairwise baselines.

The arithmetic is the core's
``sfmtool_core::geometry::translation_averaging`` kernel, reached through
``sfmtool._sfmtool.geometry.average_translations``; what is here is the
dict-keyed shape the chain speaks, converted to the kernel's arrays in sorted
key order.  See ``specs/core/geometry/translation-averaging.md``.

Provenance: the study's `relaxlib.centres_by_averaging` (316-374).  The growth
route beside it in the same module is not carried: it stalls wherever the graph
has no frontier, and the averaging route places every frame the graph connects.

What the study carried was the direction half of the objective and a linear
solve against it.  The form the objective builds sends the TRUE centres to zero
whenever the directions agree, so the constellation is the form's own null
space; the length half and the reading of that null space are what the kernel
adds.
"""

from __future__ import annotations

import numpy as np


def _kernel():
    """The core's averaging entry points, imported on use."""
    from sfmtool._sfmtool.geometry import average_translations, direction_reading

    return average_translations, direction_reading


def _arrays(frames, dirs, weights, lengths=None, length_weights=None):
    """The kernel's per-edge arrays, in sorted key order.

    ``keys`` is what the caller gets its per-edge answers back under, so the
    order every array below is built in is the one the dicts are rebuilt in."""
    index = {f: k for k, f in enumerate(frames)}
    keys = sorted(dirs)
    edges = np.array([[index[k[0]], index[k[1]]] for k in keys], np.int64).reshape(
        -1, 2
    )
    d = (
        np.ascontiguousarray(np.stack([dirs[k] for k in keys]), float)
        if keys
        else np.zeros((0, 3))
    )
    w = np.array([weights[k] for k in keys], float)
    ell = np.array([float((lengths or {}).get(k, np.nan)) for k in keys], float)
    aw = np.array([float((length_weights or {}).get(k, 0.0)) for k in keys], float)
    return keys, edges, d, w, ell, aw


def direction_reading(frames, dirs, weights):
    """What the DIRECTIONS alone determine, before any length is read in.

    The same form the averaging builds, at the weights the edges came with and
    with the length half empty, read once.  It says what the graph's geometry
    determines on its own, which is a property of the capture rather than of a
    solve, and it is what tells a colinear path from a general one."""
    _average, reading = _kernel()
    _keys, edges, d, w, _ell, _aw = _arrays(frames, dirs, weights)
    read = dict(reading(edges, d, w, len(frames)))
    # The two fields a reading does not pose: no solve was run and no length
    # was read in.
    read.pop("solved", None)
    read.pop("n_lengths", None)
    return read


def centres_by_averaging(
    frames,
    dirs,
    weights,
    lengths=None,
    length_weights=None,
    irls_rounds=None,
):
    """Camera centres from pairwise baselines, by weighted linear averaging.

    Minimizes ``sum_ij w_ij || P_ij (c_j - c_i) ||^2 + a_ij (d_ij . (c_j - c_i)
    - s L_ij)^2`` under the scale gauge ``sum_ij w_ij d_ij . (c_j - c_i) =
    sum_ij w_ij`` and the shift gauge ``sum_j c_j = 0``, with the scale ``s``
    that turns the relative lengths into distances eliminated.  Both gauges are
    exactly the freedoms the pairwise readings cannot see, so fixing them adds
    nothing.

    ``frames`` is the ordered frame list, ``dirs`` a ``{(i, j): unit d}`` in
    that frame indexing, ``weights`` a ``{(i, j): w}``.  ``lengths`` is an
    optional ``{(i, j): L}`` of relative baseline lengths on one common scale
    and ``length_weights`` how far each is trusted; an edge missing from either
    states no length and constrains only the direction.  ``irls_rounds`` of
    ``None`` takes the kernel's own count.  Returns ``(centres (n, 3),
    per-edge lambda, per-edge residual, reading)``, the first three ``None``
    where the graph states no baseline at all."""
    average, _reading = _kernel()
    keys, edges, d, w, ell, aw = _arrays(frames, dirs, weights, lengths, length_weights)
    kw = {} if irls_rounds is None else {"rounds": int(irls_rounds)}
    cen, lam, res, read = average(edges, d, w, ell, aw, n_frames=len(frames), **kw)
    read = dict(read)
    solved = bool(read.pop("solved"))
    if not solved:
        return None, None, None, read
    return cen, dict(zip(keys, lam)), dict(zip(keys, res)), read
