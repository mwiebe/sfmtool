# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 17 — can a per-point radius separate co-obs from background at all?

The cliff (and any per-point radius) can only work if, for a given descriptor,
its true co-observations are all nearer than its background neighbours. This
measures that ceiling. For every in-track descriptor that has at least one
co-observation and one background neighbour among its cross-image k-NN:

  - m_co  = distance to the *farthest* co-observation,
  - m_bg  = distance to the *nearest* background neighbour.

If ``m_co < m_bg`` the point is **separable**: a radius between them keeps all
co-obs and no background. Otherwise they interleave and no radius is clean. We
report the separability rate, how clean the gap is (``m_bg/m_co``), and the
oracle trade-off on the non-separable ones (background admitted to keep all
co-obs, vs co-obs lost to admit no background).

Usage:
    pixi run -e experiments python experiments/exp17_separability.py
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np

from exp03_radius_clusters import K
from sfm_descriptors import load_descriptor_bank

DATASETS = [
    "seoul_bull_ws",
    "seattle_backyard_ws",
    "kerry_park_ws",
    "dino_dog_toy_ws",
]


def run_one(path: str, name: str, preset: str) -> dict:
    from sfmtool import KdForest

    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    idx, dst = KdForest(desc, preset=preset).query(desc, k=K)
    idx = idx.astype(np.int64)
    img = bank.image_label.astype(np.int64)
    point = bank.point_label
    pos = point >= 0

    nb = idx[:, 1:]
    nd = dst[:, 1:]
    cross = img[nb] != img[:, None]
    co = cross & (point[nb] == point[:, None]) & pos[:, None]  # true co-obs
    bg = cross & ~((point[nb] == point[:, None]) & pos[:, None])  # cross-image bg

    has_co, has_bg = co.any(1), bg.any(1)
    src = pos & has_co & has_bg

    m_co = np.where(co, nd, -np.inf).max(1)  # farthest co-obs
    m_bg = np.where(bg, nd, np.inf).min(1)  # nearest background
    ratio = (m_bg / np.maximum(m_co, 1e-6))[src]
    separable = ratio > 1.0

    co_cnt = co[src].sum(1)
    # Oracle trade-off on the non-separable points:
    #   bg admitted if radius = m_co (keep all co-obs)
    bg_at_mco = (bg & (nd <= m_co[:, None]))[src].sum(1)
    #   co-obs lost if radius just below m_bg (admit zero background)
    co_lost_at_mbg = (co & (nd >= m_bg[:, None]))[src].sum(1)

    return dict(
        name=name,
        n_src=int(src.sum()),
        ratio=ratio,
        separable=separable,
        co_cnt=co_cnt,
        bg_at_mco=bg_at_mco,
        co_lost_at_mbg=co_lost_at_mbg,
    )


def print_table(r: dict) -> None:
    sep = r["separable"]
    print(f"\n{r['name']}: in-track points with co-obs & background = {r['n_src']}")
    print(f"  separable (a radius keeps all co-obs, no background): {sep.mean():.1%}")
    print(
        f"  median co-obs per point (≤16 NN): {np.median(r['co_cnt']):.0f}  "
        f"(p90={np.percentile(r['co_cnt'], 90):.0f})"
    )
    rs = r["ratio"][sep]
    if rs.size:
        print(
            f"  separable gap m_bg/m_co: median={np.median(rs):.2f}  "
            f"p25={np.percentile(rs, 25):.2f}"
        )
    ns = ~sep
    if ns.any():
        print(
            f"  non-separable ({ns.mean():.1%}): to keep all co-obs you admit "
            f"median {np.median(r['bg_at_mco'][ns]):.0f} background; "
            f"to admit none you lose median {np.median(r['co_lost_at_mbg'][ns]):.0f} co-obs"
        )


def plot(results, outpath):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, r in zip(axes.ravel(), results):
        lr = np.log2(np.clip(r["ratio"], 1e-3, 1e3))
        bins = np.linspace(-3, 3, 60)
        ax.hist(lr, bins=bins, color="C0", alpha=0.8)
        ax.axvline(0, color="k", ls="--", label="m_bg = m_co (boundary)")
        ax.set_title(f"{r['name']}  (separable {r['separable'].mean():.0%})")
        ax.set_xlabel("log2(m_bg / m_co)   >0: separable, <0: interleaved")
        ax.set_ylabel("count")
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Can a per-point radius separate co-obs from background?", y=0.997)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, dpi=130)
    print(f"\nwrote {outpath}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="accurate")
    ap.add_argument("--datasets", nargs="*", default=DATASETS)
    ap.add_argument(
        "--sfmr",
        default=None,
        help="explicit .sfmr file(s) to score instead of the default "
        "baseline solve (use with one or more paths)",
    )
    ap.add_argument("--sfmr-paths", nargs="*", default=None)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    print(
        "Separability ceiling: is every co-observation nearer than every "
        "background neighbour?"
    )
    explicit = args.sfmr_paths or ([args.sfmr] if args.sfmr else None)
    if explicit:
        results = [run_one(p, Path(p).stem, args.preset) for p in explicit]
    else:
        results = [
            run_one(
                sorted(glob.glob(f"../{ws}/sfmr/*solve*.sfmr"))[0],
                ws.replace("_ws", ""),
                args.preset,
            )
            for ws in args.datasets
        ]
    for r in results:
        print_table(r)
    if len(results) == len(DATASETS):
        out = Path(args.outdir)
        out.mkdir(parents=True, exist_ok=True)
        plot(results, str(out / "exp17_separability.png"))


if __name__ == "__main__":
    main()
