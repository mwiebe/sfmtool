# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0
"""Cluster strips: render a cluster-patches file's clusters for review, and
draw a blind review set of members.

A *cluster strip* is one row per cluster: the reference member's tile first,
then every measured member's tile, each sampled at the template grid of the
refinement (``resolution`` samples per side over ``patch_size`` keypoint-frame
units, from the file's ``refine_options``) at the member's refined position
and absolute shape, so every tile shows the same piece of surface as the
reference's when the member is right. Members are ordered by their scale ratio
to the reference (``sqrt|det S|`` over the reference's), and each tile is
labelled with that ratio, its status and two scores against the reference's
tile: the plain windowed ZNCC, which is what the refinement's ``min_zncc``
judges, and the blur-matched ZNCC of ``specs/core/patch/blur-matched-zncc.md``.
See ``specs/drafts/cluster-review-strips.md``.

Three subcommands::

    # Score every measured member of a file (optionally with its epipolar
    # distance under ground-truth poses), into an .npz the other two read.
    pixi run python scripts/cluster_strips.py score FILE-clusters-patches.matches \\
        -o scores.npz [--gt GROUND_TRUTH.sfmr]

    # Render chosen (or random) clusters as strips.
    pixi run python scripts/cluster_strips.py strips FILE-clusters-patches.matches \\
        --scores scores.npz (--cluster 12 --cluster 345 | --random 20) -o strips.png

    # Draw a blind review set: members rejected by the plain score and accepted
    # by the blur-matched one, stratified by scale ratio, plus controls.
    pixi run python scripts/cluster_strips.py review -o review_dir \\
        --entry NAME FILE-clusters-patches.matches scores.npz [--entry ...]

The tiles are rendered in numpy with the kernel's grid, pyramid-level rule and
pixel-centre convention, so the plain score reproduces the file's stored
``member_zncc`` (the ``score`` subcommand prints the agreement).
"""

from __future__ import annotations

import argparse
import csv
import html
import json
from pathlib import Path

import cv2
import numpy as np

from _viz_common import draw_text, new_canvas

STATUS_NAMES = {
    0: "reference",
    1: "kept",
    2: "rejected_low_zncc",
    3: "rejected_shift",
    4: "duplicate_image",
    5: "not_evaluated",
    6: "rejected_unlocalizable",
    7: "rejected_unlocalizable_refined",
    8: "rejected_unlocalizable_cells",
}
SHORT_STATUS = {0: "ref", 1: "kept", 2: "low", 3: "shift", 7: "unloc", 8: "cells"}
# Statuses whose geometry is the refinement's answer (MemberStatus::is_measured).
MEASURED = (1, 2, 3, 7, 8)
REFERENCE_UNREFINABLE = 0xFFFFFFFF
# The refinement's scoring window (ClusterRefineParams::default()); the file
# does not store it.
WINDOW_SIGMA = 0.5
# The member gates the refinement applies after the ZNCC: a rescued member is
# counted only where it would pass them too.
MAX_SHIFT_PX = 3.0
MAX_SELF_SIMILARITY_RADIUS = 2.5
# The blur-matched pair rule: the sharper tile's semi-major axis under the other
# tile's semi-minor axis, the target that semi-minor axis capped at 2 grid px,
# and at least 1.25 times the sharper tile's semi-major axis (the ratio of the
# scores against the stored bitmap).
MIN_RATIO = 1.25
MAX_TARGET = 2.0
MIN_SHARPER = 0.05
PYRAMID_LEVELS = 8
RATIO_BINS = [0.0, 0.25, 0.5, 0.7, 1.4, np.inf]
RATIO_BIN_NAMES = ["<1/4", "1/4-1/2", "1/2-0.7", "0.7-1.4", ">=1.4"]


# ===== File, workspace and images =====


def resolve_workspace(matches_path: Path) -> Path:
    """The workspace a .matches file names: its relative path, else its absolute one."""
    from sfmtool.fileio import read_matches_metadata

    meta = read_matches_metadata(str(matches_path))["workspace"]
    rel = meta.get("relative_path", "")
    if rel:
        cand = (matches_path.parent / rel).resolve()
        if (cand / ".sfm-workspace.json").exists():
            return cand
    return Path(meta["absolute_path"])


class ClusterFile:
    """A cluster-patches .matches file and the template grid it was refined at."""

    def __init__(self, path):
        from sfmtool.fileio import read_matches

        self.path = Path(path)
        d = read_matches(str(self.path))
        if not d.get("has_cluster_patches"):
            raise SystemExit(f"{path} has no cluster-patches section")
        self.workspace = resolve_workspace(self.path)
        self.names = list(d["image_names"])
        opts = d["refine_options"]
        self.radius = float(opts["patch_size"]) / 2.0
        self.resolution = int(opts["resolution"])
        self.min_zncc = float(opts["min_zncc"])
        self.starts = np.asarray(d["cluster_starts"]).astype(np.int64)
        self.member_images = np.asarray(d["member_images"]).astype(np.int64)
        self.positions = np.asarray(d["member_positions"], np.float64)
        self.shapes = np.asarray(d["member_affine_shapes"], np.float64)
        self.status = np.asarray(d["member_status"])
        self.zncc = np.asarray(d["member_zncc"])
        self.shift = np.asarray(d["member_shift_px"])
        refm = np.asarray(d["reference_members"]).astype(np.int64)
        refm[refm == REFERENCE_UNREFINABLE] = -1
        self.cluster = np.repeat(np.arange(len(self.starts) - 1), np.diff(self.starts))
        self.ref_of = refm[self.cluster]
        self.grid = TemplateGrid(self.radius, self.resolution)
        self._pyramids: dict[int, list[np.ndarray]] = {}

    def pyramid(self, image: int) -> list[np.ndarray]:
        if image not in self._pyramids:
            if len(self._pyramids) >= 16:
                self._pyramids.pop(next(iter(self._pyramids)))
            path = self.workspace / self.names[image]
            im = cv2.imread(str(path), cv2.IMREAD_COLOR)
            if im is None:
                raise SystemExit(f"cannot read {path}")
            self._pyramids[image] = box_pyramid(im)
        return self._pyramids[image]

    def drop(self, image: int) -> None:
        self._pyramids.pop(image, None)


def box_pyramid(im: np.ndarray) -> list[np.ndarray]:
    """Float levels, each the 2x2 mean of the one below (floor sizes)."""
    levels = [im.astype(np.float32)]
    while len(levels) < PYRAMID_LEVELS and min(levels[-1].shape[:2]) >= 4:
        p = levels[-1]
        h, w = p.shape[0] // 2, p.shape[1] // 2
        q = p[: 2 * h, : 2 * w]
        levels.append(
            0.25 * (q[0::2, 0::2] + q[1::2, 0::2] + q[0::2, 1::2] + q[1::2, 1::2])
        )
    return levels


# ===== The template grid =====


class TemplateGrid:
    """The refinement's R x R grid over [-radius, radius]^2 keypoint-frame units."""

    def __init__(self, radius: float, resolution: int):
        self.radius = radius
        self.res = resolution
        self.step = 2.0 * radius / resolution
        self.off = 0.5 * self.step - radius
        s = (np.arange(resolution) + 0.5) * (2.0 / resolution) - 1.0
        ss, tt = np.meshgrid(s, s)
        r2 = ss * ss + tt * tt
        self.window = np.where(
            r2 > 1.0, 0.0, np.exp(-r2 / (2 * WINDOW_SIGMA**2))
        ).ravel()
        self.support = self.window > 0

    def level(self, shapes: np.ndarray) -> np.ndarray:
        """The kernel's level rule: clamp(floor(log2 s_min), 0, L-1), with s_min
        the smaller singular value of the grid's spacing map; level 0 below 2."""
        sv = np.linalg.svd(shapes * self.step, compute_uv=False)[:, -1]
        lv = np.zeros(len(shapes), np.int64)
        ok = np.isfinite(sv) & (sv >= 2.0)
        lv[ok] = np.minimum(
            np.floor(np.log2(sv[ok])).astype(np.int64), PYRAMID_LEVELS - 1
        )
        return lv

    def render(self, pyramid, positions, shapes):
        """Tiles (n, R, R, 3) float32 and their valid flags (n, R, R): bilinear
        samples of the selected level at ``position + S u`` for each grid point
        ``u``, pixel centres at +0.5, a sample whose taps leave the frame invalid."""
        n, R = len(positions), self.res
        tiles = np.zeros((n, R, R, 3), np.float32)
        valid = np.zeros((n, R, R), bool)
        if n == 0:
            return tiles, valid
        g = self.off + self.step * np.arange(R)
        gc, gr = np.meshgrid(g, g)
        u = np.stack([gc.ravel(), gr.ravel()], 0)
        lv = self.level(shapes)
        for lev in np.unique(lv):
            sel = np.nonzero(lv == lev)[0]
            img = pyramid[min(int(lev), len(pyramid) - 1)]
            H, Wd = img.shape[:2]
            scale = float(2 ** int(lev))
            xy = positions[sel, :, None] + np.einsum("nij,jk->nik", shapes[sel], u)
            x = xy[:, 0] / scale - 0.5
            y = xy[:, 1] / scale - 0.5
            ok = (
                np.isfinite(x)
                & np.isfinite(y)
                & (x >= 0)
                & (y >= 0)
                & (x <= Wd - 1)
                & (y <= H - 1)
            )
            xc = np.clip(np.nan_to_num(x), 0, Wd - 1)
            yc = np.clip(np.nan_to_num(y), 0, H - 1)
            x0 = np.minimum(np.floor(xc).astype(np.int64), Wd - 2)
            y0 = np.minimum(np.floor(yc).astype(np.int64), H - 2)
            fx = (xc - x0)[..., None].astype(np.float32)
            fy = (yc - y0)[..., None].astype(np.float32)
            top = img[y0, x0] * (1 - fx) + img[y0, x0 + 1] * fx
            bottom = img[y0 + 1, x0] * (1 - fx) + img[y0 + 1, x0 + 1] * fx
            tiles[sel] = (top * (1 - fy) + bottom * fy).reshape(len(sel), R, R, 3)
            valid[sel] = ok.reshape(len(sel), R, R)
        return tiles, valid

    def zncc(self, a, b, va, vb):
        """Windowed ZNCC of each pair over the support samples with data in both,
        the mean over the channels with texture in ``a``; a flat ``b`` channel
        scores 0. NaN where no sample carries data in both."""
        n = len(a)
        A = a.reshape(n, -1, a.shape[-1]).astype(np.float64)
        B = b.reshape(n, -1, b.shape[-1]).astype(np.float64)
        w = self.window[None] * (
            (va.reshape(n, -1) & vb.reshape(n, -1)) & self.support[None]
        )
        sw = w.sum(1)
        wn = w[..., None] / np.maximum(sw, 1e-300)[:, None, None]
        da = A - (wn * A).sum(1, keepdims=True)
        db = B - (wn * B).sum(1, keepdims=True)
        saa = (wn * da * da).sum(1)
        sbb = (wn * db * db).sum(1)
        sab = (wn * da * db).sum(1)
        textured = saa > 1e-6
        z = np.where(
            textured & (sbb > 1e-6), sab / np.sqrt(np.maximum(saa * sbb, 1e-300)), 0.0
        )
        count = textured.sum(1)
        out = np.full(n, np.nan)
        ok = (sw > 0) & (count > 0)
        out[ok] = (z.sum(1) / np.maximum(count, 1))[ok]
        return out


# ===== Blur-matched scoring =====


def self_similarity(tiles, valid):
    """Each tile's overlap-reading self-similarity: ellipse axes (n, 2) and radius."""
    from sfmtool.patches import zncc_self_similarity_parts_stack

    if len(tiles) == 0:
        return np.zeros((0, 2)), np.zeros(0)
    stack = np.concatenate([tiles, valid[..., None].astype(np.float32) * 255.0], -1)
    r = zncc_self_similarity_parts_stack(np.ascontiguousarray(stack, np.float32))
    return np.asarray(r["ellipse_axes"]), np.asarray(r["radius"])


def pair_rule(axes_a, axes_b):
    """Which tile of each pair is blurred (0 = a, 1 = b, -1 = neither) and the
    semi-major length it is blurred to."""
    which = np.full(len(axes_a), -1)
    target = np.full(len(axes_a), np.nan)
    for side, s, o in ((0, axes_a, axes_b), (1, axes_b, axes_a)):
        major = np.maximum(s[:, 0], MIN_SHARPER)
        tgt = np.minimum(o[:, 1], MAX_TARGET)
        ok = (
            np.isfinite(major)
            & np.isfinite(tgt)
            & (s[:, 0] < o[:, 1])
            & (tgt >= MIN_RATIO * major)
            & (which < 0)
        )
        which[ok] = side
        target[ok] = tgt[ok]
    return which, target


def to_u8(t):
    return np.clip(np.rint(t), 0, 255).astype(np.uint8)


def blur_tile(tile, valid, target, assessment=None):
    """Blur one tile to the target semi-major length; returns (tile, sigma,
    assessment), the tile None where the tile's readings give no width."""
    from sfmtool.patches import assess_blur, blur_to_length

    u8 = to_u8(tile)
    if assessment is None:
        assessment = assess_blur(u8, valid=valid)
    if assessment is None:
        return None, float("nan"), None
    r = blur_to_length(u8, assessment, float(target), valid=valid)
    if r is None:
        return None, float("nan"), assessment
    return np.asarray(r["samples"], np.float32), float(r["sigma"]), assessment


def blur_matched(
    grid,
    ref_tile,
    ref_valid,
    ref_axes,
    tiles,
    valid,
    axes,
    ref_cache=None,
    ref_key=None,
):
    """The blur-matched score of each member tile against one reference tile
    each: (score, which side was blurred, sigma). ``ref_tile`` etc. are
    member-parallel."""
    plain = grid.zncc(ref_tile, tiles, ref_valid, valid)
    which, target = pair_rule(ref_axes, axes)
    score = plain.copy()
    sigma = np.full(len(tiles), np.nan)
    for k in np.nonzero(which >= 0)[0]:
        if which[k] == 0:
            key = None if ref_key is None else int(ref_key[k])
            cached = None if ref_cache is None or key is None else ref_cache.get(key)
            bt, sg, assessment = blur_tile(ref_tile[k], ref_valid[k], target[k], cached)
            if ref_cache is not None and key is not None:
                ref_cache[key] = assessment
            if bt is not None:
                score[k] = grid.zncc(
                    bt[None], tiles[k : k + 1], ref_valid[k : k + 1], valid[k : k + 1]
                )[0]
        else:
            bt, sg, _ = blur_tile(tiles[k], valid[k], target[k])
            if bt is not None:
                score[k] = grid.zncc(
                    ref_tile[k : k + 1],
                    bt[None],
                    ref_valid[k : k + 1],
                    valid[k : k + 1],
                )[0]
        sigma[k] = sg
    return plain, score, which, sigma


def scale_of(shapes):
    return np.sqrt(np.abs(np.linalg.det(shapes)))


# ===== Epipolar distance under ground-truth poses =====


def epipolar_px(gt_path, names, member_images, ref_of, positions):
    """Distance, in px of the member's photograph, of each member's position
    from the epipolar half-line of its reference's position (the directions
    from the member's camera to the points of the reference's ray in front of
    the reference's camera), under the ground truth's poses; NaN where either
    image has no pose."""
    from scipy.spatial.transform import Rotation
    from sfmtool.reconstruction import SfmrReconstruction

    rec = SfmrReconstruction.load(str(gt_path))
    gt_names = [n.replace("\\", "/").lower() for n in rec.image_names]
    gidx = np.full(len(names), -1, np.int64)
    for i, n in enumerate(names):
        n = n.replace("\\", "/").lower()
        hits = [
            j
            for j, g in enumerate(gt_names)
            if g == n or n.endswith("/" + g) or g.endswith("/" + n)
        ]
        if len(hits) == 1:
            gidx[i] = hits[0]
    q = np.asarray(rec.quaternions_wxyz)
    R = Rotation.from_quat(q[:, [1, 2, 3, 0]]).as_matrix()
    centres = -np.einsum("nij,ni->nj", R, np.asarray(rec.translations))
    cams = rec.cameras
    cam_idx = np.asarray(rec.camera_indexes)
    out = np.full(len(member_images), np.nan)
    ok = ref_of >= 0
    ok[ok] &= (gidx[member_images[ok]] >= 0) & (gidx[member_images[ref_of[ok]]] >= 0)
    idx = np.nonzero(ok)[0]
    if len(idx) == 0:
        return out
    gm = gidx[member_images[idx]]
    gr = gidx[member_images[ref_of[idx]]]

    def rays(gi, px):
        d = np.zeros((len(gi), 3))
        for ci in np.unique(cam_idx[gi]):
            s = cam_idx[gi] == ci
            d[s] = np.asarray(
                cams[ci].pixel_to_ray_batch(np.ascontiguousarray(px[s], np.float64))
            )
        d = np.einsum("nji,nj->ni", R[gi], d)
        return d / np.linalg.norm(d, axis=1, keepdims=True)

    pm = positions[idx]
    dm = rays(gm, pm)
    dr = rays(gr, positions[ref_of[idx]])
    pix = np.arccos(np.clip((dm * rays(gm, pm + [1.0, 0.0])).sum(1), -1, 1))
    w = centres[gr] - centres[gm]
    alpha = (dm * w).sum(1)
    beta = (dm * dr).sum(1)
    gamma = (w * dr).sum(1)
    ww = (w * w).sum(1)
    cands = [np.where(ww > 0, alpha / np.sqrt(np.maximum(ww, 1e-300)), -1.0), beta]
    with np.errstate(divide="ignore", invalid="ignore"):
        s_star = (alpha * gamma - beta * ww) / (beta * gamma - alpha)
        qq = ww + 2 * gamma * s_star + s_star * s_star
        cs = (alpha + beta * s_star) / np.sqrt(qq)
    cands.append(np.where(np.isfinite(s_star) & (s_star >= 0) & (qq > 0), cs, -1.0))
    out[idx] = np.arccos(np.clip(np.max(np.stack(cands), axis=0), -1, 1)) / pix
    return out


# ===== score =====


def cmd_score(args):
    cf = ClusterFile(args.matches)
    grid = cf.grid
    measured = np.isin(cf.status, MEASURED) & (cf.ref_of >= 0)
    refs = np.unique(cf.ref_of[measured])
    slot = np.full(len(cf.status), -1, np.int64)
    slot[refs] = np.arange(len(refs))
    R = grid.res
    ref_tiles = np.zeros((len(refs), R, R, 3), np.float32)
    ref_valid = np.zeros((len(refs), R, R), bool)
    for i in np.unique(cf.member_images[refs]):
        sel = refs[cf.member_images[refs] == i]
        t, v = grid.render(cf.pyramid(int(i)), cf.positions[sel], cf.shapes[sel])
        ref_tiles[slot[sel]], ref_valid[slot[sel]] = t, v
        cf.drop(int(i))
    ref_axes, _ = self_similarity(ref_tiles, ref_valid)
    M = len(cf.status)
    plain = np.full(M, np.nan)
    bm = np.full(M, np.nan)
    which = np.full(M, -2, np.int8)
    sigma = np.full(M, np.nan)
    radius = np.full(M, np.nan)
    cache: dict[int, object] = {}
    for i in range(len(cf.names)):
        sel = np.nonzero(measured & (cf.member_images == i))[0]
        if len(sel) == 0:
            continue
        t, v = grid.render(cf.pyramid(i), cf.positions[sel], cf.shapes[sel])
        cf.drop(i)
        s = slot[cf.ref_of[sel]]
        axes, radius[sel] = self_similarity(t, v)
        plain[sel], bm[sel], which[sel], sigma[sel] = blur_matched(
            grid,
            ref_tiles[s],
            ref_valid[s],
            ref_axes[s],
            t,
            v,
            axes,
            cache,
            cf.ref_of[sel],
        )
    scale = scale_of(cf.shapes)
    ratio = np.full(M, np.nan)
    ok = cf.ref_of >= 0
    ratio[ok] = scale[ok] / scale[cf.ref_of[ok]]
    epi = (
        epipolar_px(args.gt, cf.names, cf.member_images, cf.ref_of, cf.positions)
        if args.gt
        else np.full(M, np.nan)
    )
    np.savez_compressed(
        args.output,
        measured=measured,
        plain=plain,
        blur_matched=bm,
        blurred_side=which,
        blur_sigma=sigma,
        radius=radius,
        ratio=ratio,
        epipolar_px=epi,
        meta=json.dumps(
            {
                "matches": str(cf.path),
                "gt": str(args.gt) if args.gt else None,
                "window_sigma": WINDOW_SIGMA,
                "min_ratio": MIN_RATIO,
            }
        ),
    )
    d = np.abs(plain - cf.zncc)[measured]
    flips = int(((plain >= cf.min_zncc) != (cf.zncc >= cf.min_zncc))[measured].sum())
    print(
        f"{int(measured.sum())} measured members of {M}; plain score against the file's "
        f"member_zncc: |diff| median {np.nanmedian(d):.5f}, p99 {np.nanpercentile(d, 99):.5f}, "
        f"max {np.nanmax(d):.5f}, {flips} gate flips"
    )
    print(
        f"pairs blurred: reference {int((which == 0).sum())}, member {int((which == 1).sum())}"
    )


# ===== strips =====


def tile_image(tile, size):
    return cv2.resize(to_u8(tile), (size, size), interpolation=cv2.INTER_NEAREST)


def labelled_tile(tile, size, lines, colour):
    """A tile over a label plate of one text line per entry of ``lines``."""
    out = new_canvas(size, size + 14 * len(lines) + 4)
    out[:size] = tile_image(tile, size)
    cv2.rectangle(out, (0, 0), (size - 1, size - 1), colour, 1)
    for k, text in enumerate(lines):
        draw_text(out, text, (3, size + 13 + 14 * k), 0.38, (225, 225, 225))
    return out


def cluster_strip(cf, scores, c, size):
    """One cluster's strip: its reference tile, then its measured members by scale ratio."""
    members = np.arange(cf.starts[c], cf.starts[c + 1])
    ref = cf.ref_of[members[0]] if len(members) else -1
    if ref < 0:
        return None
    shown = members[np.isin(cf.status[members], MEASURED)]
    shown = shown[np.argsort(scores["ratio"][shown], kind="stable")]
    grid = cf.grid
    ids = np.concatenate([[ref], shown])
    tiles = np.zeros((len(ids), grid.res, grid.res, 3), np.float32)
    for i in np.unique(cf.member_images[ids]):
        sel = np.nonzero(cf.member_images[ids] == i)[0]
        tiles[sel], _ = grid.render(
            cf.pyramid(int(i)), cf.positions[ids[sel]], cf.shapes[ids[sel]]
        )
    blocks = [
        labelled_tile(
            tiles[0],
            size,
            ["reference", f"img {cf.member_images[ref]}", ""],
            (255, 255, 255),
        )
    ]
    for k, m in enumerate(shown, start=1):
        st = int(cf.status[m])
        p, b = scores["plain"][m], scores["blur_matched"][m]
        colour = (
            (80, 200, 80)
            if st == 1
            else ((255, 160, 60) if b >= cf.min_zncc else (70, 70, 230))
        )
        blocks.append(
            labelled_tile(
                tiles[k],
                size,
                [
                    f"x{scores['ratio'][m]:.2f} {SHORT_STATUS.get(st, str(st))}",
                    f"p {p:.3f} b {b:.3f}",
                    f"img {cf.member_images[m]}",
                ],
                colour,
            )
        )
    h = max(b.shape[0] for b in blocks)
    gap = new_canvas(4, h)
    row = []
    for b in blocks:
        pad = new_canvas(b.shape[1], h)
        pad[: b.shape[0]] = b
        row += [pad, gap]
    return np.hstack(row[:-1])


def cmd_strips(args):
    cf = ClusterFile(args.matches)
    scores = np.load(args.scores)
    clusters = list(args.cluster or [])
    if args.random:
        rng = np.random.default_rng(args.seed)
        multi = np.nonzero(
            np.bincount(cf.cluster[scores["measured"]], minlength=len(cf.starts) - 1)
            > 0
        )[0]
        clusters += rng.choice(
            multi, size=min(args.random, len(multi)), replace=False
        ).tolist()
    rows = []
    for c in clusters:
        strip = cluster_strip(cf, scores, int(c), args.tile)
        if strip is None:
            continue
        m = np.arange(cf.starts[c], cf.starts[c + 1])
        acc = m[np.isin(cf.status[m], (0, 1))]
        span = (
            scale_of(cf.shapes[acc]).max() / scale_of(cf.shapes[acc]).min()
            if len(acc) > 1
            else 1.0
        )
        label = new_canvas(150, strip.shape[0])
        for k, text in enumerate(
            [f"cluster {c}", f"{len(m)} members", f"kept span x{span:.2f}"]
        ):
            draw_text(label, text, (6, 18 + 18 * k), 0.42, (210, 210, 210))
        rows.append(np.hstack([label, new_canvas(3, strip.shape[0]), strip]))
    if not rows:
        raise SystemExit("no cluster to render")
    width = max(900, max(r.shape[1] for r in rows))
    body = []
    for r in rows:
        pad = new_canvas(width, r.shape[0])
        pad[:, : r.shape[1]] = r
        body += [pad, np.full((2, width, 3), 60, np.uint8)]
    title = new_canvas(width, 44)
    draw_text(
        title,
        f"{cf.path.name}: reference tile, then members by scale ratio (member / reference)",
        (6, 17),
        0.45,
        (255, 255, 0),
    )
    draw_text(
        title,
        "p = plain ZNCC, b = blur-matched ZNCC; frame green = kept, blue = rejected and "
        f"b >= {cf.min_zncc}, red = rejected",
        (6, 36),
        0.42,
        (255, 255, 0),
    )
    cv2.imwrite(str(args.output), np.vstack([title, *body[:-1]]))
    print(f"wrote {args.output} ({len(rows)} clusters)")


# ===== review =====

CORNERS = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)


def context_crop(cf, image, position, shape, size=288):
    """A crop of the photograph centred on the footprint, the footprint outlined
    with a tick on its top edge (the template's row direction)."""
    im = cf.pyramid(image)[0]
    quad = position + (CORNERS * cf.radius) @ shape.T
    half = max(48.0, 1.3 * float(np.abs(quad - position).max()))
    scale = size / (2 * half)
    M = np.array(
        [
            [scale, 0, size / 2 - scale * (position[0] - 0.5)],
            [0, scale, size / 2 - scale * (position[1] - 0.5)],
        ]
    )
    crop = cv2.warpAffine(
        to_u8(im),
        M,
        (size, size),
        flags=cv2.INTER_AREA if scale < 1 else cv2.INTER_CUBIC,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=(40, 40, 40),
    )

    def f(p):
        q = (p - 0.5) * scale + [
            size / 2 - scale * (position[0] - 0.5),
            size / 2 - scale * (position[1] - 0.5),
        ]
        return np.round(q * 16).astype(np.int32)

    tick = position + np.array([[0, -1.0], [0, -0.6]]) * cf.radius @ shape.T
    cv2.polylines(crop, [f(quad)], True, (0, 165, 255), 2, cv2.LINE_AA, shift=4)
    cv2.polylines(crop, [f(tick)], False, (0, 165, 255), 2, cv2.LINE_AA, shift=4)
    return crop


def case_image(cf, member, case_id, tile=288):
    """Reference and member side by side: tiles on the template grid over crops
    of their photographs with the footprints outlined. Nothing names a score or
    a status."""
    ref = int(cf.ref_of[member])
    ids = np.array([ref, member])
    tiles = []
    for m in ids:
        t, _ = cf.grid.render(
            cf.pyramid(int(cf.member_images[m])),
            cf.positions[m : m + 1],
            cf.shapes[m : m + 1],
        )
        tiles.append(tile_image(t[0], tile))
    crops = [
        context_crop(cf, int(cf.member_images[m]), cf.positions[m], cf.shapes[m], tile)
        for m in ids
    ]
    sep = new_canvas(8, tile)
    top = np.hstack([tiles[0], sep, tiles[1]])
    bottom = np.hstack([crops[0], sep, crops[1]])
    head = new_canvas(top.shape[1], 26)
    draw_text(head, f"{case_id}", (6, 18), 0.55, (255, 255, 0))
    draw_text(head, "reference", (tile // 2 - 36, 18), 0.45, (220, 220, 220))
    draw_text(head, "member", (tile + 8 + tile // 2 - 28, 18), 0.45, (220, 220, 220))
    return np.vstack([head, top, new_canvas(top.shape[1], 8), bottom])


INDEX_HTML = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cluster Member Review</title>
<style>
:root {{ --bg:#f6f6f4; --fg:#1c1c1c; --card:#fff; --line:#d8d8d4; --accent:#1f5fbf; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg:#18181a; --fg:#e8e8e6; --card:#232326; --line:#3a3a3e; --accent:#7fb0ff; }} }}
body {{ margin:0; padding:16px; background:var(--bg); color:var(--fg); font:14px/1.45 system-ui, sans-serif; }}
main {{ max-width:1200px; margin:0 auto; }}
.case {{ background:var(--card); border:1px solid var(--line); border-radius:8px; padding:12px; margin:0 0 16px; }}
.case img {{ max-width:100%; height:auto; display:block; }}
.choices label {{ margin-right:14px; }}
input[type=text] {{ width:min(600px, 100%); }}
button {{ font:inherit; padding:6px 12px; }}
#out {{ width:100%; height:120px; }}
</style></head><body><main>
<h1>Cluster member review</h1>
<p>{n} cases. Each shows a cluster's reference tile and one member's tile, both sampled on the refinement's
template grid, above crops of their photographs with the footprints outlined. For each case answer whether the
member is a <b>right</b> correspondence of the reference (the same piece of surface) or a <b>wrong</b> one.
Choices are kept in this browser while you work; <b>Export CSV</b> gives the <code>review.csv</code> rows to save
beside <code>cases.csv</code>.</p>
{cases}
<p><button id="export">Export CSV</button></p>
<textarea id="out" readonly></textarea>
</main>
<script>
const KEY = "cluster-review-" + {key};
let state = {{}};
try {{ state = JSON.parse(localStorage.getItem(KEY) || "{{}}"); }} catch (e) {{}}
function save() {{ try {{ localStorage.setItem(KEY, JSON.stringify(state)); }} catch (e) {{}} }}
document.querySelectorAll(".case").forEach(el => {{
  const id = el.dataset.id; const s = state[id] || {{}};
  el.querySelectorAll("input[type=radio]").forEach(r => {{
    if (s.choice === r.value) r.checked = true;
    r.addEventListener("change", () => {{ state[id] = Object.assign(state[id] || {{}}, {{choice: r.value}}); save(); }});
  }});
  const note = el.querySelector("input[type=text]"); note.value = s.note || "";
  note.addEventListener("input", () => {{ state[id] = Object.assign(state[id] || {{}}, {{note: note.value}}); save(); }});
}});
document.getElementById("export").addEventListener("click", () => {{
  const q = v => '"' + String(v || "").replace(/"/g, '""') + '"';
  const rows = ["case,choice,note"];
  document.querySelectorAll(".case").forEach(el => {{
    const s = state[el.dataset.id] || {{}}; rows.push([el.dataset.id, s.choice || "", q(s.note)].join(","));
  }});
  document.getElementById("out").value = rows.join("\\n");
}});
</script></body></html>
"""


def draw_review_set(entries, n_cases, n_controls, rng):
    """Draw (entry index, member, type) cases: disagreements stratified by ratio
    bin, one at a time from each (entry, bin) pool in turn; controls split
    between accepted-by-both and rejected-by-both, round robin over the
    entries."""
    pools = {}
    control_pools = {}
    for e, (cf, sc) in enumerate(entries):
        meas = sc["measured"]
        gates = (cf.shift <= MAX_SHIFT_PX) & (
            sc["radius"] <= MAX_SELF_SIMILARITY_RADIUS
        )
        low = meas & (cf.status == 2)
        accepted_bm = sc["blur_matched"] >= cf.min_zncc
        dis = low & gates & accepted_bm & (sc["blurred_side"] >= 0)
        b = np.digitize(sc["ratio"], RATIO_BINS) - 1
        for k in range(len(RATIO_BIN_NAMES)):
            pools[(e, k)] = list(rng.permutation(np.nonzero(dis & (b == k))[0]))
        control_pools[("control_accepted_both", e)] = list(
            rng.permutation(np.nonzero(meas & (cf.status == 1) & accepted_bm)[0])
        )
        control_pools[("control_rejected_both", e)] = list(
            rng.permutation(np.nonzero(low & gates & ~accepted_bm)[0])
        )

    def round_robin(pools, n, entry_of):
        out, keys = [], sorted(pools)
        while len(out) < n and any(pools[k] for k in keys):
            for k in keys:
                if pools[k] and len(out) < n:
                    out.append((entry_of(k), int(pools[k].pop())))
        return out

    picked = [
        (e, m, "disagreement") for e, m in round_robin(pools, n_cases, lambda k: k[0])
    ]
    for i, kind in enumerate(("control_accepted_both", "control_rejected_both")):
        take = n_controls // 2 + (n_controls % 2 if i == 0 else 0)
        sub = {k: v for k, v in control_pools.items() if k[0] == kind}
        picked += [(e, m, kind) for e, m in round_robin(sub, take, lambda k: k[1])]
    return [picked[i] for i in rng.permutation(len(picked))]


def cmd_review(args):
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    names = [e[0] for e in args.entry]
    entries = [(ClusterFile(m), np.load(s)) for _, m, s in args.entry]
    cases = draw_review_set(entries, args.cases, args.controls, rng)
    key_rows, blocks = [], []
    for n, (e, m, kind) in enumerate(cases, start=1):
        cf, sc = entries[e]
        cid = f"c{n:02d}"
        cv2.imwrite(str(out / f"{cid}.png"), case_image(cf, m, cid))
        ref = int(cf.ref_of[m])
        b = int(np.digitize(sc["ratio"][m], RATIO_BINS) - 1)
        key_rows.append(
            {
                "case": cid,
                "dataset": names[e],
                "type": kind,
                "member": m,
                "cluster": int(cf.cluster[m]),
                "reference_member": ref,
                "image": cf.names[cf.member_images[m]],
                "reference_image": cf.names[cf.member_images[ref]],
                "scale_ratio": f"{sc['ratio'][m]:.4f}",
                "ratio_bin": RATIO_BIN_NAMES[b],
                "status": STATUS_NAMES[int(cf.status[m])],
                "plain_zncc": f"{sc['plain'][m]:.6f}",
                "file_zncc": f"{cf.zncc[m]:.6f}",
                "blur_matched_zncc": f"{sc['blur_matched'][m]:.6f}",
                "blurred_side": {0: "reference", 1: "member"}.get(
                    int(sc["blurred_side"][m]), "none"
                ),
                "blur_sigma": f"{sc['blur_sigma'][m]:.4f}",
                "shift_px": f"{cf.shift[m]:.4f}",
                "radius": f"{sc['radius'][m]:.4f}",
                "epipolar_px": f"{sc['epipolar_px'][m]:.3f}",
                "position_x": f"{cf.positions[m, 0]:.4f}",
                "position_y": f"{cf.positions[m, 1]:.4f}",
                "s00": f"{cf.shapes[m, 0, 0]:.6f}",
                "s01": f"{cf.shapes[m, 0, 1]:.6f}",
                "s10": f"{cf.shapes[m, 1, 0]:.6f}",
                "s11": f"{cf.shapes[m, 1, 1]:.6f}",
            }
        )
        blocks.append(
            f'<section class="case" data-id="{cid}"><h2>{cid}</h2><img src="{cid}.png" alt="case {cid}">'
            f'<p class="choices"><label><input type="radio" name="{cid}" value="right"> right</label>'
            f'<label><input type="radio" name="{cid}" value="wrong"> wrong</label>'
            f'<label><input type="radio" name="{cid}" value="unsure"> unsure</label></p>'
            f'<p><input type="text" placeholder="note"></p></section>'
        )
    with open(out / "cases.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(key_rows[0]))
        w.writeheader()
        w.writerows(key_rows)
    with open(out / "review.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "choice", "note"])
        for r in key_rows:
            w.writerow([r["case"], "", ""])
    (out / "index.html").write_text(
        INDEX_HTML.format(
            n=len(key_rows),
            cases="\n".join(blocks),
            key=json.dumps(html.escape(str(out.resolve()))),
        ),
        encoding="utf-8",
    )
    kinds = {
        k: sum(r["type"] == k for r in key_rows) for k in {r["type"] for r in key_rows}
    }
    print(f"wrote {len(key_rows)} cases to {out} {kinds}; the key is cases.csv")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser(
        "score", help="score every measured member of a cluster-patches file"
    )
    p.add_argument("matches", type=Path)
    p.add_argument("-o", "--output", type=Path, required=True)
    p.add_argument(
        "--gt", type=Path, help="ground-truth .sfmr for the epipolar distance"
    )
    p = sub.add_parser("strips", help="render clusters as strips")
    p.add_argument("matches", type=Path)
    p.add_argument("--scores", type=Path, required=True)
    p.add_argument("--cluster", type=int, action="append")
    p.add_argument("--random", type=int, default=0)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--tile", type=int, default=100)
    p.add_argument("-o", "--output", type=Path, required=True)
    p = sub.add_parser("review", help="draw a blind review set")
    p.add_argument(
        "--entry",
        nargs=3,
        action="append",
        required=True,
        metavar=("NAME", "MATCHES", "SCORES"),
    )
    p.add_argument("--cases", type=int, default=40)
    p.add_argument("--controls", type=int, default=8)
    p.add_argument("--seed", type=int, default=20261010)
    p.add_argument("-o", "--output", type=Path, required=True)
    args = ap.parse_args()
    {"score": cmd_score, "strips": cmd_strips, "review": cmd_review}[args.cmd](args)


if __name__ == "__main__":
    main()
