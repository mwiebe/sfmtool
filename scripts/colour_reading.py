# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0
"""Colour reading: the pairwise colour agreement of a cluster member's tile
with its reference's, after the smallest photographic difference between the
two is fitted out.

Step 1 of ``specs/drafts/image-photometric-model.md``. For each member it
renders the reference's tile ``A`` and the member's tile ``B`` on the
refinement's template grid exactly as ``cluster_strips.py`` does, blur-matches
the pair by the pair rule of ``specs/core/patch/blur-matched-zncc.md`` (the
sharper tile is blurred), and then:

1. **Linearises** both tiles with one gamma, ``(v / 255) ** gamma``, for each
   gamma of ``--gamma`` (no response curve is known; both entries are one
   camera each, so one gamma serves both tiles), and divides both by ``A``'s
   mean intensity over the samples it fits, so a residual is in units of
   ``A``'s normalised radiance.
2. **Masks** clipped samples. A photograph pixel is clipped when any channel is
   0 or 255. Clipping is read on each photograph's own pixels: the clipped
   pixel masks are reduced with the same box pyramid as the photograph and
   sampled through the same grid, so a sample is clipped when any photograph
   pixel under its bilinear taps is; and the tile's clipped share is the share
   of full-resolution pixel centres inside its outline (the grid's square at
   the member's position and shape) that are clipped, the glossary's
   **clipped share**.
3. **Fits** ``B -> A`` over the unmasked samples of the scoring window's
   support, weighted by the window: an exposure difference ``dt`` (``B`` is
   scaled by ``2 ** dt``), PPISP's chromaticity homography on RG chromaticity
   and intensity built from four offsets ``dc_k`` of the R, G, B and W
   primaries (``H = T diag(k) S^-1``, ``H33 = 1``, intensity preserved), and a
   desaturation ``s`` in ``[0, 0.99]`` toward luminance. The cost is the
   window-weighted mean over samples and channels of a Huber loss on
   ``A - fit(B)``, plus the priors at ``--prior-weight`` (PPISP's weight is
   1): a Huber loss with delta 0.1 on ``dt`` and 0.005 on the length of each
   ``dc_k``. PPISP puts these terms on the mean over all frames, the gauge;
   here they act on the pair's difference, where at weight 1 they outweigh
   the data term (a 1.5x exposure fits as ``2 ** 0.44``), so weight 0 is the
   other setting measured. ``dt`` is bounded to +-8 stops. The data
   term's delta is read from the pair: a least-squares fit first, then the
   Huber fit from it with delta 1.345 times the MAD scale of its residuals.
   **Rank rule**: the chromaticity scatter of each tile's unmasked samples
   (window-weighted standard deviations along its principal axes) is read, and
   the smaller tile's is kept; both singular values at least ``RANK_TOL`` fits
   the full model (10 parameters), only the first fits a white-point shift
   (``dc_W`` alone, which is a per-channel gain with intensity kept, plus
   ``dt``), neither fits ``dt`` alone.
4. **Reads** the colour residual (the window-weighted RMS of the per-sample
   residual over the unmasked samples, the worst 10% by weight trimmed, in
   ``A``'s normalised radiance; the untrimmed RMS too), a chromaticity residual
   (the same over the samples' RG chromaticities), and the plausibility of the
   fit: ``|dt|`` in stops, the largest ``|dc_k|``, the white-point disagreement
   (the RMS distance of the four fitted target chromaticities from the best
   per-channel-gain homography's, which moves only W), ``s``, and the clip
   coincidence (of the samples clipped in exactly one tile, the share whose
   other tile is in the matching extreme quartile of its intensities: top for a
   255, bottom for a 0). The tone of each tile, over every sample with data
   including clipped ones: B's brightness over A's in stops and each tile's
   saturation, so a glare-like pair (one tile brighter, less saturated and
   clipped) can be counted whichever tile it is on. Beside them: the plain and blur-matched ZNCC of the
   pair, and a baseline, the residual of a per-channel gain fit (3 gains,
   Huber, no priors).

Two subcommands::

    pixi run python scripts/colour_reading.py measure FILE-clusters-patches.matches \\
        SCORES.npz --name ENTRY -o readings.npz [--include cases.csv] \\
        [--accepted-per-bin 2000 --accepted-min 5000] [--gamma 1.0 --gamma 2.2]

    pixi run python scripts/colour_reading.py review -o review_dir \\
        --entry NAME MATCHES SCORES READINGS [--entry ...] --cases 20

``measure`` reads the ``score`` subcommand's ``.npz`` of ``cluster_strips.py``
for the populations and the epipolar distances; ``review`` draws the members
on which the colour residual and the blur-matched ZNCC disagree, rendered in
the cluster-strips convention (``cluster_strips.case_image``), with
``cases.csv`` the key.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path

import cv2
import numpy as np
from scipy.optimize import minimize

from cluster_strips import (
    MAX_SELF_SIMILARITY_RADIUS,
    MAX_SHIFT_PX,
    RATIO_BIN_NAMES,
    RATIO_BINS,
    STATUS_NAMES,
    ClusterFile,
    blur_tile,
    box_pyramid,
    case_image,
    pair_rule,
    self_similarity,
)

HUBER_K = 1.345
DT_DELTA = 0.1
DC_DELTA = 0.005
RANK_TOL = 0.01
TRIM = 0.10
MIN_SAMPLES = 20
S_MAX = 0.99
DT_MAX = 8.0
LUMA = np.array([0.2126, 0.7152, 0.0722])
# RGB -> (R, G, I = R + G + B), PPISP's RGI, and back.
C = np.array([[1.0, 0, 0], [0, 1.0, 0], [1.0, 1.0, 1.0]])
CINV = np.linalg.inv(C)
# Homogeneous source chromaticities of R, G, B (columns) and W.
SRC = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 1.0]])
SRC_W = np.array([1 / 3, 1 / 3, 1.0])
SRC_PTS = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0], [1 / 3, 1 / 3]])


# ===== The model =====


def homography(dc):
    """PPISP's chromaticity homography from offsets ``dc`` (4, 2) of R, G, B, W."""
    t = SRC_PTS + dc
    T = np.vstack([t[:3].T, np.ones(3)])
    k = np.linalg.solve(T, np.array([t[3, 0], t[3, 1], 1.0]))
    H = T @ np.diag(k) @ np.linalg.inv(SRC @ np.diag(np.linalg.solve(SRC, SRC_W)))
    return H / H[2, 2]


def apply_model(B, dt, H, s):
    """``fit(B)``: exposure, chromaticity homography with intensity kept, desaturation."""
    x = (2.0**dt) * B @ C.T
    y = x @ H.T
    n = x[:, 2] / np.where(np.abs(y[:, 2]) > 1e-9, y[:, 2], 1e-9)
    rgb = (y * n[:, None]) @ CINV.T
    if s:
        rgb = (1 - s) * rgb + s * (rgb @ LUMA)[:, None]
    return rgb


def huber(x, d):
    a = np.abs(x)
    return np.where(a <= d, 0.5 * x * x, d * (a - 0.5 * d))


def unpack(p, rank):
    dt = p[0]
    dc = np.zeros((4, 2))
    s = 0.0
    if rank == 2:
        dc = p[1:9].reshape(4, 2)
        s = p[9]
    elif rank == 1:
        dc[3] = p[1:3]
    return dt, dc, s


def huber_scale(r):
    """The Huber delta for residuals ``r``: 1.345 times their MAD scale."""
    return max(HUBER_K * 1.4826 * float(np.median(np.abs(r))), 1e-4)


def fit_pair(A, B, w, rank, prior_weight=1.0):
    """Fit ``B -> A``; returns (dt, dc (4, 2), s, delta). A least-squares fit
    first, then a Huber fit from it with the delta read off its residuals."""
    wn = w / w.sum()
    ia, ib = (wn * A.sum(1)).sum(), (wn * B.sum(1)).sum()
    npar = {0: 1, 1: 3, 2: 10}[rank]
    p = np.zeros(npar)
    p[0] = np.log2(max(ia, 1e-6) / max(ib, 1e-6))
    bounds = [(-DT_MAX, DT_MAX)] + [(None, None)] * (npar - 1)
    if rank == 2:
        bounds[9] = (0.0, S_MAX)

    def cost(q, delta):
        dt, dc, s = unpack(q, rank)
        try:
            H = homography(dc)
        except np.linalg.LinAlgError:
            return 1e6
        r = A - apply_model(B, dt, H, s)
        c = (wn[:, None] * huber(r, delta)).sum() / 3.0
        c += prior_weight * huber(dt, DT_DELTA)
        c += prior_weight * huber(np.linalg.norm(dc, axis=1), DC_DELTA).sum()
        return c if np.isfinite(c) else 1e6

    p = minimize(cost, p, args=(np.inf,), method="L-BFGS-B", bounds=bounds).x
    dt, dc, s = unpack(p, rank)
    delta = huber_scale(A - apply_model(B, dt, homography(dc), s))
    p = minimize(cost, p, args=(delta,), method="L-BFGS-B", bounds=bounds).x
    return (*unpack(p, rank), delta)


def gain_fit(A, B, w, iters=5):
    """Per-channel gains: least squares, then Huber by IRLS with the delta read
    off the least-squares residuals; the fitted ``B``."""
    ww = np.repeat(w[:, None], 3, 1)
    g = (ww * A * B).sum(0) / np.maximum((ww * B * B).sum(0), 1e-12)
    delta = huber_scale(A - B * g)
    for _ in range(iters):
        r = np.abs(A - B * g)
        ww = w[:, None] * np.where(r <= delta, 1.0, delta / np.maximum(r, 1e-12))
        g = (ww * A * B).sum(0) / np.maximum((ww * B * B).sum(0), 1e-12)
    return B * g


def residual(A, F, w):
    """(trimmed RMS, RMS) of the per-sample residual, window-weighted."""
    e = ((A - F) ** 2).mean(1)
    order = np.argsort(e)
    cw = np.cumsum(w[order]) / w.sum()
    keep = order[cw <= 1 - TRIM + 1e-12]
    if len(keep) == 0:
        keep = order[:1]
    trimmed = np.sqrt((w[keep] * e[keep]).sum() / w[keep].sum())
    return float(trimmed), float(np.sqrt((w * e).sum() / w.sum()))


def chroma(X):
    i = np.maximum(X.sum(1), 1e-9)
    return np.stack([X[:, 0] / i, X[:, 1] / i], 1)


def chroma_scatter(X, w):
    """Window-weighted standard deviations of the RG chromaticities along their
    principal axes, largest first."""
    c = chroma(X)
    wn = w / w.sum()
    d = c - (wn[:, None] * c).sum(0)
    cov = (wn[:, None, None] * d[:, :, None] * d[:, None, :]).sum(0)
    return np.sqrt(np.maximum(np.linalg.eigvalsh(cov)[::-1], 0.0))


def white_point_disagreement(dc):
    """RMS distance, over the four primaries, of the fitted targets from those
    of the best per-channel-gain homography. A gain homography keeps R, G and B
    where they are and moves only W, so the best one puts W on the fitted W and
    the distance is that of the R, G and B offsets."""
    return float(np.sqrt((dc[:3] ** 2).sum() / 4.0))


def read_pair(a8, b8, va, vb, clip_a, clip_b, window, support, gamma, prior_weight=1.0):
    """All readings of one blur-matched pair at one gamma; a dict."""
    out = {}
    m = va & vb & support & ~clip_a.any(1) & ~clip_b.any(1)
    out["n_samples"] = int(m.sum())
    if m.sum() < MIN_SAMPLES:
        return out
    w = window[m]
    A = (np.clip(a8[m], 0, 255) / 255.0) ** gamma
    B = (np.clip(b8[m], 0, 255) / 255.0) ** gamma
    ia = (w * A.mean(1)).sum() / w.sum()
    if ia <= 1e-6:
        return out
    A, B = A / ia, B / ia
    sa, sb = chroma_scatter(A, w), chroma_scatter(B, w)
    sv = np.minimum(sa, sb)
    rank = 2 if sv[1] >= RANK_TOL else (1 if sv[0] >= RANK_TOL else 0)
    t0 = time.perf_counter()
    dt, dc, s, delta = fit_pair(A, B, w, rank, prior_weight)
    out["huber_delta"] = delta
    out["fit_seconds"] = time.perf_counter() - t0
    F = apply_model(B, dt, homography(dc), s)
    out["resid"], out["resid_rms"] = residual(A, F, w)
    out["chroma_resid"], _ = residual(chroma(A), chroma(F), w)
    out["gain_resid"], out["gain_resid_rms"] = residual(A, gain_fit(A, B, w), w)
    out["rank"] = rank
    out["sv0"], out["sv1"] = float(sv[0]), float(sv[1])
    out["dt"] = float(dt)
    out["dc_max"] = float(np.linalg.norm(dc, axis=1).max())
    out["dc_w"] = float(np.linalg.norm(dc[3]))
    out["wp_disagreement"] = white_point_disagreement(dc)
    out["s"] = float(s)
    out["dc"] = dc.ravel().tolist()
    return out


def tile_tone(a8, b8, va, vb, window, support, gamma):
    """Brightness and saturation of each tile over every sample with data in
    both, clipped ones included: (log2 of B's mean intensity over A's, A's
    saturation, B's saturation), the saturation the window-weighted mean
    distance of the RG chromaticity from white, (1/3, 1/3)."""
    m = va & vb & support
    if m.sum() < MIN_SAMPLES:
        return float("nan"), float("nan"), float("nan")
    w = window[m]
    A = (np.clip(a8[m], 0, 255) / 255.0) ** gamma
    B = (np.clip(b8[m], 0, 255) / 255.0) ** gamma
    ia, ib = (w * A.sum(1)).sum(), (w * B.sum(1)).sum()
    sat = [
        float((w * np.linalg.norm(chroma(X) - 1 / 3, axis=1)).sum() / w.sum())
        for X in (A, B)
    ]
    return float(np.log2(max(ib, 1e-9) / max(ia, 1e-9))), sat[0], sat[1]


def clip_coincidence(a8, b8, va, vb, clip_a, clip_b, support):
    """Of the samples clipped in exactly one tile, the share whose other tile is
    in the matching extreme quartile of its own intensities (top for a 255,
    bottom for a 0); NaN where none is."""
    m = va & vb & support
    ia, ib = a8.mean(1), b8.mean(1)
    hits = total = 0
    for clip, other, oclip in ((clip_b, ia, clip_a), (clip_a, ib, clip_b)):
        ref = m & ~oclip.any(1)
        if ref.sum() < 4:
            continue
        lo, hi = np.percentile(other[ref], [25, 75])
        only = m & clip.any(1) & ~oclip.any(1)
        hi_c = only & clip[:, 0]
        lo_c = only & ~clip[:, 0] & clip[:, 1]
        hits += int((hi_c & (other >= hi)).sum() + (lo_c & (other <= lo)).sum())
        total += int(hi_c.sum() + lo_c.sum())
    return hits / total if total else float("nan")


# ===== Clipping on the photograph =====


class ClipImage:
    """One photograph's clipped-pixel masks: a box pyramid of (high, low, 0)
    shares for the grid's sampler, and a row-wise prefix count of clipped pixels
    for the clipped share inside an outline."""

    def __init__(self, im):
        high = (im == 255).any(2)
        low = (im == 0).any(2)
        stack = np.stack([high, low, np.zeros_like(high)], 2).astype(np.float32)
        self.pyramid = box_pyramid(stack)
        clipped = (high | low).astype(np.int32)
        self.prefix = np.concatenate(
            [np.zeros((im.shape[0], 1), np.int32), np.cumsum(clipped, 1)], 1
        )
        self.h, self.w = im.shape[:2]

    def share(self, outline):
        """Share of pixel centres inside a convex outline that are clipped."""
        y0 = max(int(np.floor(outline[:, 1].min())), 0)
        y1 = min(int(np.ceil(outline[:, 1].max())), self.h - 1)
        if y1 < y0:
            return float("nan")
        yc = np.arange(y0, y1 + 1) + 0.5
        p, q = outline, np.roll(outline, -1, 0)
        lo = np.full(len(yc), np.inf)
        hi = np.full(len(yc), -np.inf)
        for (xi, yi), (xj, yj) in zip(p, q):
            cross = (yi > yc) != (yj > yc)
            x = np.where(cross, xi + (xj - xi) * (yc - yi) / ((yj - yi) or 1e-12), 0)
            lo = np.where(cross, np.minimum(lo, x), lo)
            hi = np.where(cross, np.maximum(hi, x), hi)
        ok = np.isfinite(lo)
        a = np.clip(np.ceil(lo[ok] - 0.5), 0, self.w).astype(np.int64)
        b = np.clip(np.ceil(hi[ok] - 0.5), 0, self.w).astype(np.int64)
        b = np.maximum(a, b)
        rows = np.arange(y0, y1 + 1)[ok]
        total = int((b - a).sum())
        if total == 0:
            cx = np.clip(outline[:, 0].astype(np.int64), 0, self.w - 1)
            cy = np.clip(outline[:, 1].astype(np.int64), 0, self.h - 1)
            v = self.prefix[cy, cx + 1] - self.prefix[cy, cx]
            return float(v.mean())
        clipped = int((self.prefix[rows, b] - self.prefix[rows, a]).sum())
        return clipped / total


CORNERS = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)


# ===== measure =====


def select_members(
    cf, sc, rng, per_bin, minimum, include, only_include=False, rejected=False
):
    """Rescued members (all), a ratio-stratified sample of accepted members,
    and the members named by ``include``; (members, population) arrays."""
    meas = sc["measured"]
    gates = (cf.shift <= MAX_SHIFT_PX) & (sc["radius"] <= MAX_SELF_SIMILARITY_RADIUS)
    low = meas & (cf.status == 2)
    rescued = (
        low & gates & (sc["blur_matched"] >= cf.min_zncc) & (sc["blurred_side"] >= 0)
    )
    kept = meas & (cf.status == 1)
    if rejected:
        # Rejected low-ZNCC members past the later gates whose blur-matched
        # score also misses the bar, in place of the accepted sample.
        kept = low & gates & (sc["blur_matched"] < cf.min_zncc)
        rescued = np.zeros_like(rescued)
    b = np.digitize(sc["ratio"], RATIO_BINS) - 1
    picks = []
    left = []
    for k in range(len(RATIO_BIN_NAMES)):
        pool = rng.permutation(np.nonzero(kept & (b == k))[0])
        picks.append(pool[:per_bin])
        left.append(pool[per_bin:])
    n = sum(len(p) for p in picks)
    if n < minimum:
        big = int(np.argmax([len(x) for x in left]))
        picks.append(left[big][: minimum - n])
    acc = np.concatenate(picks)
    members = {int(m): "rescued" for m in np.nonzero(rescued)[0]}
    for m in acc:
        members.setdefault(int(m), "rejected" if rejected else "accepted")
    if only_include:
        members = {m: v for m, v in members.items() if m in set(include)}
    for m in include:
        members.setdefault(int(m), "review")
    ids = np.array(sorted(members), np.int64)
    pop = np.array([members[int(m)] for m in ids])
    return ids, pop


def cmd_measure(args):
    cf = ClusterFile(args.matches)
    sc = np.load(args.scores)
    grid = cf.grid
    rng = np.random.default_rng(args.seed)
    include = []
    if args.include:
        with open(args.include, newline="") as f:
            include = [
                int(r["member"]) for r in csv.DictReader(f) if r["dataset"] == args.name
            ]
    ids, pop = select_members(
        cf,
        sc,
        rng,
        args.accepted_per_bin,
        args.accepted_min,
        include,
        args.only_include,
        args.rejected,
    )
    refs = np.unique(cf.ref_of[ids])
    print(
        f"{args.name}: {len(ids)} members ({dict(zip(*np.unique(pop, return_counts=True)))}), {len(refs)} references"
    )

    def render_all(members):
        R = grid.res
        tiles = np.zeros((len(members), R, R, 3), np.float32)
        valid = np.zeros((len(members), R, R), bool)
        clip = np.zeros((len(members), R, R, 2), bool)
        share = np.full(len(members), np.nan)
        imgs = cf.member_images[members]
        for i in np.unique(imgs):
            sel = np.nonzero(imgs == i)[0]
            m = members[sel]
            tiles[sel], valid[sel] = grid.render(
                cf.pyramid(int(i)), cf.positions[m], cf.shapes[m]
            )
            ci = ClipImage(cf.pyramid(int(i))[0].astype(np.uint8))
            ct, _ = grid.render(ci.pyramid, cf.positions[m], cf.shapes[m])
            clip[sel] = ct[..., :2] > 1e-6
            for k, mm in zip(sel, m):
                share[k] = ci.share(
                    cf.positions[mm] + (CORNERS * cf.radius) @ cf.shapes[mm].T
                )
            cf.drop(int(i))
        return tiles, valid, clip, share

    t0 = time.perf_counter()
    rt, rv, rc, rshare = render_all(refs)
    mt, mv, mc, mshare = render_all(ids)
    slot = np.searchsorted(refs, cf.ref_of[ids])
    ref_axes, _ = self_similarity(rt, rv)
    mem_axes, _ = self_similarity(mt, mv)
    which, target = pair_rule(ref_axes[slot], mem_axes)
    print(f"rendered in {time.perf_counter() - t0:.1f} s")
    gammas = list(args.gamma)
    R2 = grid.res * grid.res
    rows: list[dict] = []
    cache: dict[int, object] = {}
    t_start = time.perf_counter()
    for k in range(len(ids)):
        r = int(slot[k])
        a, b = rt[r], mt[k]
        t_pair = time.perf_counter()
        if which[k] == 0:
            bt, _, ass = blur_tile(a, rv[r], target[k], cache.get(r))
            cache[r] = ass
            if bt is not None:
                a = bt
        elif which[k] == 1:
            bt, _, _ = blur_tile(b, mv[k], target[k])
            if bt is not None:
                b = bt
        t_blur = time.perf_counter() - t_pair
        a8, b8 = a.reshape(R2, 3), b.reshape(R2, 3)
        va, vb = rv[r].ravel(), mv[k].ravel()
        ca, cb = rc[r].reshape(R2, 2), mc[k].reshape(R2, 2)
        sup = grid.support
        row = {
            "clip_share_a": rshare[r],
            "clip_share_b": mshare[k],
            "clip_samples_a": float(
                (ca.any(1) & va & sup).sum() / max(1, (va & sup).sum())
            ),
            "clip_samples_b": float(
                (cb.any(1) & vb & sup).sum() / max(1, (vb & sup).sum())
            ),
            "clip_coincidence": clip_coincidence(a8, b8, va, vb, ca, cb, sup),
            "zncc_plain_here": float(
                grid.zncc(rt[r][None], mt[k][None], rv[r][None], mv[k][None])[0]
            ),
            "zncc_bm_here": float(
                grid.zncc(a[None], b[None], rv[r][None], mv[k][None])[0]
            ),
            "blur_seconds": t_blur,
        }
        for g in gammas:
            (
                row[f"g{g:g}_tone_stops"],
                row[f"g{g:g}_sat_a"],
                row[f"g{g:g}_sat_b"],
            ) = tile_tone(a8, b8, va, vb, grid.window, sup, g)
            if args.no_fit:
                continue
            rd = read_pair(
                a8, b8, va, vb, ca, cb, grid.window, sup, g, args.prior_weight
            )
            for key, v in rd.items():
                row[f"g{g:g}_{key}"] = v
        row["pair_seconds"] = time.perf_counter() - t_pair
        rows.append(row)
        if (k + 1) % 1000 == 0:
            el = time.perf_counter() - t_start
            print(f"  {k + 1}/{len(ids)} pairs, {el:.0f} s", flush=True)
    keys = sorted({key for row in rows for key in row if not key.endswith("_dc")})
    cols = {key: [row.get(key, np.nan) for row in rows] for key in keys}
    dcs = {
        f"g{g:g}_dc": np.array(
            [row.get(f"g{g:g}_dc", [np.nan] * 8) for row in rows], np.float64
        )
        for g in gammas
    }
    out = {key: np.asarray(v, np.float64) for key, v in cols.items()}
    np.savez_compressed(
        args.output,
        member=ids,
        population=pop,
        status=cf.status[ids],
        cluster=cf.cluster[ids],
        reference=cf.ref_of[ids],
        ratio=sc["ratio"][ids],
        plain=sc["plain"][ids],
        blur_matched=sc["blur_matched"][ids],
        blurred_side=which,
        epipolar_px=sc["epipolar_px"][ids],
        meta=json.dumps(
            {
                "matches": str(cf.path),
                "scores": str(args.scores),
                "name": args.name,
                "gammas": gammas,
                "huber_k": HUBER_K,
                "dt_delta": DT_DELTA,
                "dc_delta": DC_DELTA,
                "rank_tol": RANK_TOL,
                "trim": TRIM,
                "seed": args.seed,
                "prior_weight": args.prior_weight,
                "dt_max": DT_MAX,
            }
        ),
        **out,
        **dcs,
    )
    print(
        f"wrote {args.output}: {len(ids)} pairs in {time.perf_counter() - t_start:.0f} s"
    )


# ===== review =====


def cmd_review(args):
    """Members where the colour residual and the blur-matched ZNCC disagree:
    the residual above the median of the entry's epipolar-wrong accepted and
    rescued members while the blur-matched ZNCC reaches 0.85, and the residual
    at most the median of the epipolar-right ones while the blur-matched ZNCC
    misses 0.85. ``READINGS`` may name several ``measure`` outputs, comma
    separated (a ``--rejected`` run supplies the second kind). Up to ``--cases`` of each kind, round robin over the entries, then
    shuffled and rendered in the cluster-strips convention."""
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    key = f"g{args.gamma_key:g}_resid"
    pools = {}
    entries = []
    for e, (name, matches, scores, readings) in enumerate(args.entry):
        cf = ClusterFile(matches)
        sc = np.load(scores)
        parts = [np.load(r) for r in readings.split(",")]
        rd = {
            k: np.concatenate([q[k] for q in parts])
            for k in parts[0].files
            if k != "meta" and all(k in q.files for q in parts)
        }
        entries.append((name, cf, sc, rd))
        epi = rd["epipolar_px"]
        kept = rd["status"] == 1
        bar = max(
            3.0,
            5
            * float(np.nanmedian(sc["epipolar_px"][sc["measured"] & (cf.status == 1)])),
        )
        wrong = np.isfinite(epi) & (epi > bar)
        right = np.isfinite(epi) & (epi <= bar)
        res = rd[key]
        hi = float(
            np.nanmedian(
                res[wrong & np.isin(rd["population"], ["accepted", "rescued"])]
            )
        )
        lo = float(
            np.nanmedian(
                res[right & np.isin(rd["population"], ["accepted", "rescued"])]
            )
        )
        bm = rd["blur_matched"]
        sel = np.isfinite(res) & (rd["population"] != "review")
        pools[("colour_wrong_zncc_right", e)] = list(
            rng.permutation(rd["member"][sel & (res > hi) & (bm >= 0.85)])
        )
        pools[("colour_right_zncc_wrong", e)] = list(
            rng.permutation(rd["member"][sel & (res <= lo) & (bm < 0.85)])
        )
        print(
            f"{name}: residual bars wrong-median {hi:.4f}, right-median {lo:.4f}; pools "
            f"{len(pools[('colour_wrong_zncc_right', e)])} / {len(pools[('colour_right_zncc_wrong', e)])}"
            f" (kept {int(kept.sum())})"
        )
    picked = []
    for kind in ("colour_wrong_zncc_right", "colour_right_zncc_wrong"):
        sub = [k for k in sorted(pools) if k[0] == kind]
        got = 0
        while got < args.cases and any(pools[k] for k in sub):
            for k in sub:
                if pools[k] and got < args.cases:
                    picked.append((k[1], int(pools[k].pop()), kind))
                    got += 1
    picked = [picked[i] for i in rng.permutation(len(picked))]
    rows = []
    for n, (e, m, kind) in enumerate(picked, start=1):
        name, cf, sc, rd = entries[e]
        cid = f"d{n:02d}"
        cv2.imwrite(str(out / f"{cid}.png"), case_image(cf, m, cid))
        j = int(np.nonzero(rd["member"] == m)[0][0])
        rows.append(
            {
                "case": cid,
                "dataset": name,
                "type": kind,
                "member": m,
                "cluster": int(cf.cluster[m]),
                "reference_member": int(cf.ref_of[m]),
                "image": cf.names[cf.member_images[m]],
                "reference_image": cf.names[cf.member_images[cf.ref_of[m]]],
                "scale_ratio": f"{sc['ratio'][m]:.4f}",
                "ratio_bin": RATIO_BIN_NAMES[
                    int(np.digitize(sc["ratio"][m], RATIO_BINS) - 1)
                ],
                "status": STATUS_NAMES[int(cf.status[m])],
                "population": str(rd["population"][j]),
                "plain_zncc": f"{sc['plain'][m]:.6f}",
                "blur_matched_zncc": f"{sc['blur_matched'][m]:.6f}",
                "colour_residual": f"{rd[key][j]:.5f}",
                "gain_residual": f"{rd[f'g{args.gamma_key:g}_gain_resid'][j]:.5f}",
                "epipolar_px": f"{sc['epipolar_px'][m]:.3f}",
            }
        )
    with open(out / "cases.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(out / "review.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "choice", "note"])
        for r in rows:
            w.writerow([r["case"], "", ""])
    kinds = {k: sum(r["type"] == k for r in rows) for k in {r["type"] for r in rows}}
    print(f"wrote {len(rows)} cases to {out} {kinds}; the key is cases.csv")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("measure", help="read the colour of chosen members")
    p.add_argument("matches", type=Path)
    p.add_argument("scores", type=Path)
    p.add_argument("--name", required=True)
    p.add_argument("--include", type=Path, help="a cases.csv whose members to add")
    p.add_argument(
        "--no-fit",
        action="store_true",
        help="the tone and clipping readings only, no model fit",
    )
    p.add_argument(
        "--rejected",
        action="store_true",
        help="sample rejected low-ZNCC members the blur-matched score also rejects",
    )
    p.add_argument(
        "--only-include", action="store_true", help="read the --include members only"
    )
    p.add_argument("--accepted-per-bin", type=int, default=2000)
    p.add_argument("--accepted-min", type=int, default=5000)
    p.add_argument("--gamma", type=float, action="append")
    p.add_argument(
        "--prior-weight",
        type=float,
        default=1.0,
        help="weight of the priors on dt and dc_k (PPISP's is 1)",
    )
    p.add_argument("--seed", type=int, default=20261010)
    p.add_argument("-o", "--output", type=Path, required=True)
    p = sub.add_parser(
        "review", help="draw a blind review set where colour and ZNCC disagree"
    )
    p.add_argument(
        "--entry",
        nargs=4,
        action="append",
        required=True,
        metavar=("NAME", "MATCHES", "SCORES", "READINGS"),
    )
    p.add_argument("--cases", type=int, default=20)
    p.add_argument("--gamma-key", type=float, default=2.2)
    p.add_argument("--seed", type=int, default=20261011)
    p.add_argument("-o", "--output", type=Path, required=True)
    args = ap.parse_args()
    if args.cmd == "measure" and not args.gamma:
        args.gamma = [1.0, 2.2]
    {"measure": cmd_measure, "review": cmd_review}[args.cmd](args)


if __name__ == "__main__":
    main()
