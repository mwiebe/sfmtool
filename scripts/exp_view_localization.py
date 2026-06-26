#!/usr/bin/env python3
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment: how well-localized is each view's projection onto a surfel?

The good-view experiments grew the refinement view set and raised the consensus
Phi, but the normal's *determinacy* never sharpened. One suspect is residual
**mis-registration**: refinement assumes a view's projection lands exactly on the
patch (the warp comes straight from the point's 3D position + pose). If the point
center is a little off, or a pose/distortion residual shifts the projection, each
view's rendered patch is translated by a fraction of a pixel to a few pixels
relative to the consensus -- which blurs the consensus template and flattens the
photo-consistency peak even when the *normal* is right.

This probe measures that directly. For each view of a well-supported point:

1. Render it over a **larger context tile** than the scored core (so the search
   window can slide without running off the patch edge -- the boundary problem).
2. Build a **leave-one-out** robust consensus reference from the *other* views
   (so a view is never aligned to a template it polluted).
3. **Coarse-to-fine, sub-pixel translation search** (downsample for the coarse
   pass, parabolic refine at full res) for the shift that maximizes windowed
   ZNCC to that reference.

It then reports, per point: the per-view shift magnitude (in patch pixels and in
*source* pixels), how much the ZNCC gains from aligning, the consensus Phi
centered vs. per-view-aligned (the localization headroom), and a decomposition of
the shifts into a **common** in-plane component (the point center is off -- a
re-centering would fix every view at once) vs. the **residual** incoherent spread
(per-view pose/distortion/depth error the search can't blame on the center). A
re-render at the found offset (translating the patch center, a single resample)
verifies the gain is real and not an artifact of resampling the context tile.

Usage::

    pixi run python scripts/exp_view_localization.py RECON.sfmr -o out.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import cv2
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _normal_strip_lib import (  # noqa: E402
    ANG_RANGE,
    EXTENT_FACTOR,
    INIT_STEPS,
    PATCH,
    _phi,
    gauss_window,
    geometric_views,
    irls_weights,
    legend_bar,
    render_cores,
    vet,
    znorm_stack,
)

from sfmtool._solve_strips import _SolveStrips  # noqa: E402
from sfmtool._sfmtool import OrientedPatch, PatchCloud, SfmrReconstruction  # noqa: E402
from sfmtool._sfmtool.flow import WarpMap  # noqa: E402

TRACK_COL = (0, 255, 255)  # yellow: feature-matched (track) view
GOOD_COL = (90, 230, 90)  # green: admitted by photometric agreement


def _znorm_img(img: np.ndarray, w: np.ndarray) -> np.ndarray:
    """z-normalized (C, P) vector for one patch image, matching the lib's
    windowed convention (sqrt(w) folded in so plain dots are windowed ZNCCs)."""
    return znorm_stack([img], w)[0]


def _zncc(a: np.ndarray, b: np.ndarray) -> float:
    """Windowed ZNCC of two z-normalized (C, P) patches, mean over channels."""
    return float((a * b).sum() / a.shape[0])


def _parabolic(mid: float, left: float, right: float) -> float:
    """Sub-sample peak offset in [-0.5, 0.5] from a 3-point quadratic fit
    (scores at -1, 0, +1 around an integer maximum)."""
    denom = left - 2.0 * mid + right
    if abs(denom) < 1e-12:
        return 0.0
    return float(np.clip(0.5 * (left - right) / denom, -1.0, 1.0))


def _window(tile: np.ndarray, off: int, dy: int, dx: int) -> np.ndarray:
    return tile[off + dy : off + dy + PATCH, off + dx : off + dx + PATCH]


def _subpix_window(tile: np.ndarray, off: int, dy: float, dx: float) -> np.ndarray:
    """A PATCH-sized window centered at the (fractional) offset, bilinearly."""
    cx = off + (PATCH - 1) / 2.0 + dx
    cy = off + (PATCH - 1) / 2.0 + dy
    return cv2.getRectSubPix(tile.astype(np.float32), (PATCH, PATCH), (cx, cy))


def search_shift(tile, ref_img, off, w, search, coarse):
    """Coarse-to-fine sub-pixel translation search of ``tile``'s core window
    against ``ref_img``. Returns ``(dx, dy, c_center, c_peak)`` -- the best shift
    in patch pixels and the ZNCC at the centered vs. found position.

    The coarse pass runs on a half-res pyramid level (the "downsampling": it
    covers the full ``search`` range at a quarter of the candidate count and is
    robust to local texture noise); the fine pass is a full-res integer refine
    plus a separable parabolic sub-pixel fit."""
    lo, hi = -search, search
    # bound the window to the rendered context so it never reads past the edge
    lo = max(lo, -off)
    hi = min(hi, tile.shape[0] - PATCH - off)
    base = _znorm_img(ref_img, w)

    seed_y = seed_x = 0
    if coarse and (hi - lo) >= 4:
        t2 = cv2.pyrDown(tile.astype(np.float32))
        r2 = cv2.pyrDown(ref_img.astype(np.float32))
        p2 = PATCH // 2
        off2 = (t2.shape[0] - p2) // 2
        w2 = gauss_window(p2)
        b2 = _znorm_img(r2, w2)
        best = (-2.0, 0, 0)
        rc = (hi - lo) // 2 // 2  # half-res steps spanning the full range
        for dy in range(-rc, rc + 1):
            for dx in range(-rc, rc + 1):
                y, x = off2 + dy, off2 + dx
                if y < 0 or x < 0 or y + p2 > t2.shape[0] or x + p2 > t2.shape[1]:
                    continue
                c = _zncc(_znorm_img(t2[y : y + p2, x : x + p2], w2), b2)
                if c > best[0]:
                    best = (c, dy * 2, dx * 2)
        seed_y, seed_x = best[1], best[2]

    # Full-res integer refine in a small window around the coarse seed.
    fine = 2 if coarse else search
    surf = {}
    for dy in range(seed_y - fine, seed_y + fine + 1):
        for dx in range(seed_x - fine, seed_x + fine + 1):
            if not (lo <= dy <= hi and lo <= dx <= hi):
                continue
            surf[(dy, dx)] = _zncc(_znorm_img(_window(tile, off, dy, dx), w), base)
    (py, px) = max(surf, key=surf.get)
    c_peak = surf[(py, px)]

    # Separable parabolic sub-pixel using the integer neighbors (fall back to the
    # integer peak at a search-window edge where a neighbor is missing).
    def nb(dy, dx):
        return surf.get((dy, dx))

    sy = sx = 0.0
    if nb(py - 1, px) is not None and nb(py + 1, px) is not None:
        sy = _parabolic(c_peak, nb(py - 1, px), nb(py + 1, px))
    if nb(py, px - 1) is not None and nb(py, px + 1) is not None:
        sx = _parabolic(c_peak, nb(py, px - 1), nb(py, px + 1))

    c_center = surf.get((0, 0))
    if c_center is None:  # centered offset fell outside the refine window
        c_center = _zncc(_znorm_img(_window(tile, off, 0, 0), w), base)
    return px + sx, py + sy, c_center, c_peak


def loo_reference(cores, rows, v):
    """Leave-one-out robust consensus *image* (H, W, C) from every view but ``v``
    -- the reference view ``v`` is aligned to (it never sees its own pixels)."""
    idx = [k for k in range(len(cores)) if k != v]
    raw = np.stack([np.asarray(cores[k], np.float64) for k in idx])
    wt = (
        irls_weights(rows[idx])
        if len(idx) >= 2
        else np.ones(len(idx)) / max(len(idx), 1)
    )
    return (wt[:, None, None, None] * raw).sum(0)


def _patch_of(s, pid, normal, views, ext):
    center = s.positions[int(pid)]
    up = s.rot_of[views[0]].T @ np.array([0.0, -1.0, 0.0])
    return center, OrientedPatch.from_center_normal(
        center.tolist(), np.asarray(normal).tolist(), up.tolist(), [ext, ext]
    )


def patch_frame(s, pid, normal, views, ext, render_res):
    """Patch center, *unit* in-plane axes, and world-units-per-patch-pixel -- the
    frame needed to translate the patch by a shift expressed in patch pixels."""
    center, patch = _patch_of(s, pid, normal, views, ext)
    u_ax = np.asarray(patch.u_axis, float)
    u_ax /= np.linalg.norm(u_ax)
    v_ax = np.asarray(patch.v_axis, float)
    v_ax /= np.linalg.norm(v_ax)
    return center, u_ax, v_ax, 2.0 * ext / render_res


def render_tile_at_offset(
    s, i, center, normal, u_ax, v_ax, dx, dy, wpp, ext, render_res
):
    """Full context tile for view ``i`` with the patch center translated by
    ``(dx, dy)`` patch px in-plane -- a *single* resample of the source image (no
    re-sampling of an already-warped tile, so applying a shift this way never
    compounds blur)."""
    shifted = center + (dx * wpp) * u_ax + (dy * wpp) * v_ax
    up = s.rot_of[i].T @ np.array([0.0, -1.0, 0.0])
    patch = OrientedPatch.from_center_normal(
        shifted.tolist(), np.asarray(normal).tolist(), up.tolist(), [ext, ext]
    )
    wm = WarpMap.from_patch(patch, s.cam_of[i], s.pose_of[i], render_res)
    return np.asarray(wm.remap_bilinear(s.image(i)), np.float32)


def select_good_view_points(
    recon, s, imgs, render_res, w, common, *, num, pool, accept_frac, floor
):
    """Pick ``num`` well-supported points and grow each a good-view set + normal.

    The shared front end for the localization probes: candidate points are the
    finite, well-tracked ones with extra geometrically-visible views; each gets a
    track-only seed normal, a photometrically-vetted good-view set, and a normal
    re-refined over that set. Returns dicts ``{pid, good, track, normal, ext}``
    ranked by how many extra views they admit."""
    finite = np.abs(np.asarray(recon.positions_xyzw)[:, 3]) > 1e-9
    cidx = {int(p): k for k, p in enumerate(s.cloud.point_ids)}

    cand = []
    for pid in s.cloud.point_ids:
        pid = int(pid)
        if not finite[pid]:
            continue
        tracked = sorted(set(s.obs.get(pid, [])))
        if len(tracked) < 3:
            continue
        vis = geometric_views(s, pid)
        if set(vis) - set(tracked):
            cand.append((pid, set(tracked), vis, len(set(vis) - set(tracked))))
    if not cand:
        raise SystemExit("no finite, well-tracked points with extra visible views")
    cand.sort(key=lambda t: t[3], reverse=True)
    pool_c = cand[:pool]
    pool_ids = [pid for pid, *_ in pool_c]
    ext_of = {pid: s._half(pid) * (render_res / PATCH) for pid in pool_ids}

    cloud_t = PatchCloud.from_reconstruction(
        recon, normal="stored", extent_value=EXTENT_FACTOR
    )
    res_t = cloud_t.refine_normals(recon, imgs, point_ids=pool_ids, **common)
    n_track = {pid: np.asarray(res_t["normal"])[cidx[pid]] for pid in pool_ids}

    scored = []
    for pid, T, vis, _ in pool_c:
        _, _, _, admitted, _ = vet(
            s,
            pid,
            n_track[pid],
            vis,
            T,
            sorted(T),
            w,
            render_res,
            ext_of[pid],
            accept_frac,
            floor,
        )
        scored.append((pid, T, sorted(T | set(admitted)), len(admitted)))
    scored.sort(key=lambda r: r[3], reverse=True)
    sel = scored[:num]
    sel_ids = [r[0] for r in sel]

    view_idx = [[] for _ in s.cloud.point_ids]
    for pid in sel_ids:
        view_idx[cidx[pid]] = [int(i) for i in dict(((r[0], r[2]) for r in sel))[pid]]
    cloud_g = PatchCloud.from_reconstruction(
        recon, normal="stored", extent_value=EXTENT_FACTOR
    )
    res_g = cloud_g.refine_normals(
        recon, imgs, point_ids=sel_ids, view_indices=view_idx, **common
    )

    return [
        {
            "pid": pid,
            "good": good,
            "track": T,
            "ext": ext_of[pid],
            "normal": np.asarray(res_g["normal"])[cidx[pid]],
        }
        for pid, T, good, _ in sel
    ]


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("recon", type=Path, help="solved .sfmr")
    ap.add_argument("-o", "--out", type=Path, default=Path("view_localization.png"))
    ap.add_argument("--num", type=int, default=8, help="points to analyze")
    ap.add_argument(
        "--montage", type=int, default=4, help="points drawn in the montage"
    )
    ap.add_argument(
        "--context",
        type=int,
        default=64,
        help="context render size (>= PATCH); "
        "the search slides the PATCH core inside this without hitting the edge",
    )
    ap.add_argument(
        "--search", type=int, default=10, help="max shift searched (patch px)"
    )
    ap.add_argument(
        "--no-coarse", action="store_true", help="full-res search only (no pyramid)"
    )
    ap.add_argument("--pool", type=int, default=120, help="candidate points pre-vetted")
    ap.add_argument("--accept-frac", type=float, default=0.7)
    ap.add_argument("--floor", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    recon = SfmrReconstruction.load(str(args.recon))
    s = _SolveStrips(
        recon, recon.workspace_dir, patch=PATCH, extent_factor=EXTENT_FACTOR
    )
    imgs = [s.image(i) for i in range(len(s.names))]
    render_res = max(args.context, PATCH)
    off = (render_res - PATCH) // 2
    w = gauss_window(PATCH)
    common = dict(resolution=PATCH, angular_range_deg=ANG_RANGE, init_steps=INIT_STEPS)
    coarse = not args.no_coarse

    pts = select_good_view_points(
        recon,
        s,
        imgs,
        render_res,
        w,
        common,
        num=args.num,
        pool=args.pool,
        accept_frac=args.accept_frac,
        floor=args.floor,
    )
    sel_ids = [p["pid"] for p in pts]
    good_of = {p["pid"]: p["good"] for p in pts}
    track_of = {p["pid"]: p["track"] for p in pts}
    n_good = {p["pid"]: p["normal"] for p in pts}
    ext_of = {p["pid"]: p["ext"] for p in pts}

    print(
        f"\n{'pt':>6} {'views':>5} {'shift(src)':>11} {'ZNCC gain':>10} "
        f"{'Phi c->aln':>16} {'common/resid':>13} {'exact':>7}"
    )
    print("-" * 76)

    montage_rows = []
    agg = {"shift": [], "gain": [], "dphi": [], "common": [], "resid": [], "exact": []}
    for pid in sel_ids:
        G = good_of[pid]
        normal = n_good[pid]
        ext = ext_of[pid]
        tiles, cores = render_cores(s, pid, normal, G, render_res, ext)
        rows = znorm_stack(cores, w)

        # World scale + in-plane axes for converting patch-px shifts to source px.
        center, patch = _patch_of(s, pid, normal, G, ext)
        u_ax = np.asarray(patch.u_axis, float)
        u_ax /= np.linalg.norm(u_ax)
        v_ax = np.asarray(patch.v_axis, float)
        v_ax /= np.linalg.norm(v_ax)
        wpp = 2.0 * ext / render_res  # world units per patch pixel

        shifts, gains, src_px, aligned_rows = [], [], [], []
        exact_ok = 0
        for v, i in enumerate(G):
            ref = loo_reference(cores, rows, v)
            dx, dy, c0, c1 = search_shift(tiles[v], ref, off, w, args.search, coarse)
            shifts.append((dx, dy))
            gains.append(c1 - c0)
            # source-pixel displacement of the found shift in this view
            shifted = center + (dx * wpp) * u_ax + (dy * wpp) * v_ax
            src_px.append(_proj_dist(s, i, center, shifted))
            # aligned row for the per-view-aligned consensus
            aw = _subpix_window(tiles[v], off, dy, dx)
            aligned_rows.append(_znorm_img(aw, w))
            # exact verification: re-render at the offset (one resample) and compare
            ec = _exact_zncc(
                s,
                i,
                center,
                normal,
                u_ax,
                v_ax,
                dx,
                dy,
                wpp,
                ext,
                render_res,
                off,
                w,
                _znorm_img(ref, w),
            )
            exact_ok += int(ec >= c1 - 0.03)  # exact confirms (within slack) the gain

        aligned = np.asarray(aligned_rows)
        wt0 = irls_weights(rows)
        phi0 = _phi(rows, wt0)
        wt1 = irls_weights(aligned)
        phi1 = _phi(aligned, wt1)

        sv = np.asarray(shifts)  # (V, 2) in (dx, dy) patch px
        common_vec = sv.mean(0)
        resid = sv - common_vec
        common_px = float(np.hypot(*common_vec)) * wpp_to_src(
            s, G, center, u_ax, v_ax, wpp
        )
        resid_px = float(np.sqrt((resid**2).sum(1).mean())) * wpp_to_src(
            s, G, center, u_ax, v_ax, wpp
        )
        med_src = float(np.median(src_px))
        med_gain = float(np.median(gains))

        agg["shift"].append(med_src)
        agg["gain"].append(med_gain)
        agg["dphi"].append(phi1 - phi0)
        agg["common"].append(common_px)
        agg["resid"].append(resid_px)
        agg["exact"].append(exact_ok / len(G))

        print(
            f"{pid:>6} {len(G):>5} {med_src:>8.2f}px  {med_gain:>+9.3f}  "
            f"{phi0:>+6.3f}->{phi1:<+6.3f}  {common_px:>4.2f}/{resid_px:<4.2f}  "
            f"{exact_ok:>2}/{len(G)}"
        )

        if len(montage_rows) < args.montage:
            montage_rows.append(
                _montage_row(tiles, G, track_of[pid], shifts, gains, off, pid)
            )

    print("-" * 86)
    print(
        f"median over {len(sel_ids)} pts: shift {np.median(agg['shift']):.2f}px src, "
        f"ZNCC gain {np.median(agg['gain']):+.3f}, dPhi {np.median(agg['dphi']):+.3f}, "
        f"common {np.median(agg['common']):.2f}px vs resid {np.median(agg['resid']):.2f}px, "
        f"exact-confirm {100 * np.mean(agg['exact']):.0f}%"
    )

    if montage_rows:
        width = max(r.shape[1] for r in montage_rows)
        montage_rows = [_pad_to(r, width) for r in montage_rows]
        legend = legend_bar(
            width,
            [
                (TRACK_COL, "track view"),
                (GOOD_COL, "good view"),
                ((255, 255, 255), "core extent"),
                ((60, 200, 255), "found shift"),
            ],
        )
        out = np.vstack([legend, *montage_rows])
        args.out.parent.mkdir(parents=True, exist_ok=True)
        cv2.imwrite(str(args.out), out)
        print(f"\nwrote {args.out}")


def wpp_to_src(s, views, center, u_ax, v_ax, wpp) -> float:
    """Mean source px per 1 patch-px of in-plane shift, over ``views`` -- so a
    shift in patch px can be reported in interpretable source pixels."""
    d = []
    for i in views:
        d.append(_proj_dist(s, i, center, center + wpp * u_ax))
    return float(np.mean(d)) if d else 1.0


def _proj_dist(s, i, c0, c1) -> float:
    def proj(c):
        pc = s.rot_of[i] @ (c - s.centers[i])
        if pc[2] <= 1e-6:
            return None
        return np.array(s.cam_of[i].project(pc[0] / pc[2], pc[1] / pc[2]))

    a, b = proj(c0), proj(c1)
    if a is None or b is None:
        return float("nan")
    return float(np.hypot(*(b - a)))


def _exact_zncc(
    s, i, center, normal, u_ax, v_ax, dx, dy, wpp, ext, render_res, off, w, ref_unit
):
    """Re-render the patch with its center translated by the found in-plane shift
    (a single resample of the source) and ZNCC the core to the reference -- if the
    context-tile search were fooled by double resampling, this would disagree."""
    shifted = center + (dx * wpp) * u_ax + (dy * wpp) * v_ax
    up = s.rot_of[i].T @ np.array([0.0, -1.0, 0.0])
    patch = OrientedPatch.from_center_normal(
        shifted.tolist(), np.asarray(normal).tolist(), up.tolist(), [ext, ext]
    )
    wm = WarpMap.from_patch(patch, s.cam_of[i], s.pose_of[i], render_res)
    full = np.asarray(wm.remap_bilinear(s.image(i)), np.float32)
    core = full[off : off + PATCH, off : off + PATCH]
    return _zncc(_znorm_img(core, w), ref_unit)


def _montage_row(tiles, views, track, shifts, gains, off, pid, disp=120, gutter=150):
    sep = np.full((disp, 2, 3), 40, np.uint8)
    cells = []
    render_res = tiles[0].shape[0]
    scale = disp / render_res
    for t, i, (dx, dy), g in zip(tiles, views, shifts, gains):
        bgr = cv2.resize(
            np.clip(t, 0, 255).astype(np.uint8),
            (disp, disp),
            interpolation=cv2.INTER_NEAREST,
        )
        col = TRACK_COL if i in track else GOOD_COL
        # core extent box (white) and the found, shifted window (cyan)
        x0, x1 = round(off * scale), round((off + PATCH) * scale)
        cv2.rectangle(bgr, (x0, x0), (x1 - 1, x1 - 1), (255, 255, 255), 1)
        sx0 = round((off + dx) * scale)
        sy0 = round((off + dy) * scale)
        sx1 = round((off + dx + PATCH) * scale)
        sy1 = round((off + dy + PATCH) * scale)
        cv2.rectangle(bgr, (sx0, sy0), (sx1 - 1, sy1 - 1), (60, 200, 255), 1)
        cv2.rectangle(bgr, (0, 0), (disp - 1, disp - 1), col, 3)
        _chip(bgr, str(i), (4, 14), col)
        _chip(
            bgr,
            f"{np.hypot(dx, dy):.1f}px {g:+.2f}",
            (4, disp - 7),
            (60, 200, 255),
            0.32,
        )
        cells.extend((bgr, sep))
    strip = np.hstack(cells[:-1]) if cells else np.zeros((disp, disp, 3), np.uint8)
    g = np.full((disp, gutter, 3), 25, np.uint8)
    cv2.putText(
        g,
        f"pt {pid}",
        (6, disp // 2),
        cv2.FONT_HERSHEY_SIMPLEX,
        0.5,
        (230, 230, 230),
        1,
        cv2.LINE_AA,
    )
    return np.hstack([g, strip])


def _chip(img, text, org, color, scale=0.34) -> None:
    (tw, th), _ = cv2.getTextSize(text, cv2.FONT_HERSHEY_SIMPLEX, scale, 1)
    x, y = org
    cv2.rectangle(img, (x - 1, y - th - 2), (x + tw + 1, y + 2), (15, 15, 15), -1)
    cv2.putText(
        img, text, (x, y), cv2.FONT_HERSHEY_SIMPLEX, scale, color, 1, cv2.LINE_AA
    )


def _pad_to(row, width):
    if row.shape[1] >= width:
        return row
    pad = np.full((row.shape[0], width - row.shape[1], 3), 25, np.uint8)
    return np.hstack([row, pad])


if __name__ == "__main__":
    main()
