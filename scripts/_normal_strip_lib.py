# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Shared helpers for the normal-refinement strip experiments.

Patch rendering, the spec's z-normalized robust consensus
(``specs/core/patch-normal-refinement.md``), and strip drawing — used by both
``exp_mvs_normal_strips.py`` (refine over *all* visible views) and
``exp_goodview_normal_strips.py`` (refine over the track, then grow the view set
by photometric agreement). These run the real Rust refiner via
``PatchCloud.refine_normals(view_indices=...)``; the per-view weights and
correlations drawn on the strips are computed here, in Python, on exactly the
tiles displayed.
"""

from __future__ import annotations

import cv2
import numpy as np

from sfmtool._sfmtool import OrientedPatch
from sfmtool._sfmtool.flow import WarpMap

PATCH = 32  # validated surfel extent (px), as in compare --strips
EXTENT_FACTOR = 5.0  # PatchCloud extent as a multiple of feature size
ANG_RANGE = 25.0  # normal-search angular half-range (deg)
INIT_STEPS = 7  # coarse-grid steps per axis


def geometric_views(s, pid: int) -> list[int]:
    """Every image index whose camera geometrically sees point ``pid``: the point
    projects in front of the camera and its centre lands inside the frame. This is
    the MVS view set — a superset of the point's feature track (it ignores
    self-occlusion, which the photometric tests downstream are there to catch)."""
    c = s.positions[int(pid)]
    out: list[int] = []
    for i in range(len(s.names)):
        pc = s.rot_of[i] @ (c - s.centers[i])
        if pc[2] <= 1e-6:  # behind the camera
            continue
        cam = s.cam_of[i]
        u, v = cam.project(pc[0] / pc[2], pc[1] / pc[2])
        if 0.0 <= u < cam.width and 0.0 <= v < cam.height:
            out.append(i)
    return out


def gauss_window(n: int) -> np.ndarray:
    u = np.arange(n) - n / 2 + 0.5
    gx, gy = np.meshgrid(u, u)
    return np.exp(-(gx**2 + gy**2) / (2 * (n / 4.0) ** 2)).ravel()


def znorm_stack(cores: list[np.ndarray], w: np.ndarray) -> np.ndarray:
    """Per-view, per-channel z-normalized patch vectors with the scoring window
    folded in. Each colour channel is centred and unit-normalized under the window
    inner product ``<a,b>_w = Σ w·a·b`` (the spec's z-normalization), then scaled
    by ``sqrt(w)`` so plain Euclidean dot products equal the windowed ones.
    Returns an ``(V, C, P)`` array (views × channels × pixels), each ``[v, c]``
    row a unit vector."""
    g = np.sqrt(w)
    out = []
    for p in cores:
        a = p.reshape(-1, p.shape[-1]) if p.ndim == 3 else p.reshape(-1, 1)
        chans = []
        for c in range(a.shape[1]):
            x = a[:, c].astype(np.float64)
            x = x - (w * x).sum() / w.sum()
            nrm = np.sqrt((w * x * x).sum())
            chans.append(g * (x / nrm if nrm > 1e-9 else np.zeros_like(x)))
        out.append(np.stack(chans, 0))
    return np.asarray(out)


def _template(rows: np.ndarray, wt: np.ndarray) -> np.ndarray:
    """Weighted consensus mean ``x̄_w = Σ wᵢ xᵢ`` — shape ``(C, P)``."""
    return (wt[:, None, None] * rows).sum(0)


def _phi(rows: np.ndarray, wt: np.ndarray) -> float:
    """Weighted mean-pairwise consensus ρ̄_w, averaged over channels."""
    xbar = _template(rows, wt)
    s2 = float((wt * wt).sum())
    denom = 1.0 - s2
    if denom <= 1e-9:
        return float("nan")
    nbar2 = float((xbar**2).sum(1).mean())  # mean over channels of ‖x̄^(c)‖²
    return (nbar2 - s2) / denom


def irls_weights(rows: np.ndarray, iters: int = 5) -> np.ndarray:
    """Robust per-view inclusion weights (spec ``RobustWeighted``): Tukey weights
    on each view's residual ``‖xᵢ − x̄_w‖`` with a MAD scale, re-formed ``iters``
    times. Sums to 1."""
    v = rows.shape[0]
    if v < 2:
        return np.ones(v)
    wt = np.full(v, 1.0 / v)
    for _ in range(iters):
        xbar = _template(rows, wt)
        res = np.sqrt(((rows - xbar) ** 2).sum((1, 2)))
        scale = 1.4826 * np.median(res) if v > 2 else res.mean() + 1e-9
        c = 4.685 * scale if scale > 1e-9 else 1e-9
        u = np.clip(res / c, 0, 1)
        tk = (1 - u * u) ** 2
        if tk.sum() < 1e-9:
            tk = np.ones(v)
        wt = tk / tk.sum()
    return wt


def irls_consensus(cores: list[np.ndarray], w: np.ndarray, iters: int = 5):
    """``(weights, phi)`` over a patch stack — convenience over the steps above."""
    rows = znorm_stack(cores, w)
    if rows.shape[0] < 2:
        return np.ones(rows.shape[0]), float("nan")
    wt = irls_weights(rows, iters)
    return wt, _phi(rows, wt)


def corr_to_template(rows: np.ndarray, subset: list[int]):
    """Per-view correlation of every view to the (robust) consensus template built
    from ``subset`` (indices into ``rows``). Returns ``(corr, phi_subset)`` where
    ``corr[v]`` is the mean-over-channels ZNCC of view ``v`` against the subset's
    unit template, and ``phi_subset`` is the subset's own consensus."""
    sub = znorm_stack_subset(rows, subset)
    wt = irls_weights(sub) if len(subset) >= 2 else np.ones(len(subset))
    xbar = _template(sub, wt)
    nrm = np.sqrt((xbar**2).sum(1, keepdims=True))
    m = xbar / np.where(nrm > 1e-9, nrm, 1.0)  # unit template per channel, (C, P)
    corr = (rows * m[None]).sum((1, 2)) / rows.shape[1]  # mean over channels
    return corr, _phi(sub, wt)


def znorm_stack_subset(rows: np.ndarray, subset: list[int]) -> np.ndarray:
    return rows[np.asarray(subset, dtype=int)]


def vet(s, pid, normal, A, T, good, w, render_res, ext, accept_frac, floor):
    """Render point ``pid`` under ``normal`` over views ``A``, correlate each to
    the consensus template built from ``good`` (⊇ track ``T``), and decide which
    extra views to admit. Returns ``(rows, corr, thr, admitted, base)`` where
    ``corr`` is the per-view ZNCC to the template, ``thr`` the adaptive admit
    threshold, and ``base`` the track's own mean self-correlation."""
    _, cores = render_cores(s, pid, normal, A, render_res, ext)
    rows = znorm_stack(cores, w)
    good_pos = [A.index(i) for i in good]
    corr, _ = corr_to_template(rows, good_pos)
    track_pos = [A.index(i) for i in sorted(T)]
    base = float(np.mean(corr[track_pos]))
    thr = max(floor, accept_frac * base)
    admitted = [A[k] for k in range(len(A)) if A[k] not in T and corr[k] >= thr]
    return rows, corr, thr, admitted, base


def render_cores(s, pid: int, normal, views, render_res: int, ext: float):
    """Render point ``pid`` into each view under ``normal`` at ``render_res``;
    returns the full padded tiles and the inner ``PATCH``-sized validated cores."""
    center = s.positions[int(pid)]
    up = s.rot_of[views[0]].T @ np.array([0.0, -1.0, 0.0])
    patch = OrientedPatch.from_center_normal(
        center.tolist(), np.asarray(normal).tolist(), up.tolist(), [ext, ext]
    )
    off = (render_res - PATCH) // 2
    tiles, cores = [], []
    for i in views:
        wm = WarpMap.from_patch(patch, s.cam_of[i], s.pose_of[i], render_res)
        full = np.asarray(wm.remap_bilinear(s.image(i)), np.float32)
        tiles.append(full)
        cores.append(full[off : off + PATCH, off : off + PATCH] if off > 0 else full)
    return tiles, cores


def _chip_text(img, text, org, color, scale=0.34) -> None:
    """Draw ``text`` on a small dark chip so it stays legible over busy texture."""
    (tw, th), _ = cv2.getTextSize(text, cv2.FONT_HERSHEY_SIMPLEX, scale, 1)
    x, y = org
    cv2.rectangle(img, (x - 1, y - th - 2), (x + tw + 1, y + 2), (15, 15, 15), -1)
    cv2.putText(img, text, (x, y), cv2.FONT_HERSHEY_SIMPLEX, scale, color, 1, cv2.LINE_AA)


def strip_image(
    tiles,
    views,
    idx_colors,
    notes,
    bar_vals,
    render_res: int,
    *,
    label: str,
    disp: int = 110,
    gutter: int = 150,
) -> np.ndarray:
    """One horizontal strip of padded tiles. Each tile gets a thick border in its
    category colour ``idx_colors[i]`` (the primary included/excluded signal), a
    white box around the validated extent, its image index on a dark chip, and a
    bottom quality bar + ``notes[i]`` (green→red by ``bar_vals[i]``). A tile whose
    ``bar_vals[i]`` is ``None`` is *excluded*: it is dimmed and labelled, not
    scored — so excluded views visibly recede behind a coloured frame."""
    off = (render_res - PATCH) // 2
    sep = np.full((disp, 2, 3), 40, np.uint8)
    row: list[np.ndarray] = []
    for t, i, col, note, bar in zip(tiles, views, idx_colors, notes, bar_vals):
        p8 = np.clip(t, 0, 255).astype(np.uint8)
        bgr = p8 if p8.ndim == 3 else cv2.cvtColor(p8, cv2.COLOR_GRAY2BGR)
        scale = disp / bgr.shape[0]
        bgr = cv2.resize(bgr, (disp, disp), interpolation=cv2.INTER_NEAREST)
        excluded = bar is None
        if excluded:
            bgr = (bgr * 0.4).astype(np.uint8)  # dim views that didn't make the cut
        if off > 0:
            x0, x1 = round(off * scale), round((off + PATCH) * scale)
            cv2.rectangle(bgr, (x0, x0), (x1 - 1, x1 - 1), (255, 255, 255), 1)
        # Thick category border — the at-a-glance included/excluded marker.
        cv2.rectangle(bgr, (0, 0), (disp - 1, disp - 1), col, 4)
        _chip_text(bgr, str(i), (4, 15), col)
        if bar is not None:
            c = (int(60 * (1 - bar)), int(220 * bar + 30), int(220 * (1 - bar)))
            cv2.rectangle(bgr, (5, disp - 9), (5 + int(bar * (disp - 10)), disp - 5), c, -1)
            _chip_text(bgr, note, (4, disp - 12), c, scale=0.32)
        else:
            _chip_text(bgr, "excluded", (4, disp - 8), col, scale=0.32)
        row.extend((bgr, sep))
    strip = np.hstack(row[:-1]) if row else np.zeros((disp, disp, 3), np.uint8)
    g = np.full((disp, gutter, 3), 25, np.uint8)
    cv2.putText(g, label, (6, disp // 2), cv2.FONT_HERSHEY_SIMPLEX, 0.42, (230, 230, 230), 1, cv2.LINE_AA)
    return np.hstack([g, strip])


def legend_bar(width: int, items, *, height: int = 26) -> np.ndarray:
    """A legend row: coloured swatches + labels for ``items`` = ``[(bgr, text)]``."""
    bar = np.full((height, width, 3), 30, np.uint8)
    x = 8
    for color, text in items:
        cv2.rectangle(bar, (x, 6, 18, height - 12), color, -1)
        cv2.rectangle(bar, (x, 6, 18, height - 12), (200, 200, 200), 1)
        x += 24
        (tw, _), _ = cv2.getTextSize(text, cv2.FONT_HERSHEY_SIMPLEX, 0.42, 1)
        cv2.putText(bar, text, (x, height - 8), cv2.FONT_HERSHEY_SIMPLEX, 0.42, (225, 225, 225), 1, cv2.LINE_AA)
        x += tw + 26
    return bar
