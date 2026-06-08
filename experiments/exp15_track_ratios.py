# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 15 — neighbour-distance ratios for real tracks.

To choose the cliff radius rule `keep neighbours within c·d2` (equivalently
`dK/d1 ≥ c·(d2/d1)` marks the first excluded), look at where true co-observations
actually sit. For every in-track descriptor we label its k-NN neighbours as
same-track cross-image co-observations or background, and measure each neighbour's
distance as a multiple of `d2` (the 2nd-nearest). A good `c` keeps most
co-observations while admitting little background.

Reports, per dataset: the `d2/d1` distribution; co-observation and background
distance/`d2` distributions; and a sweep of `c` showing co-observation recall vs
background admitted. Also writes a histogram.

Usage:
    pixi run -e experiments python experiments/exp15_track_ratios.py
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
CS = [1.25, 1.5, 2.0, 2.5, 3.0, 4.0]


def run_one(ws: str, preset: str) -> dict:
    from sfmtool import KdForest

    path = sorted(glob.glob(f"../{ws}/sfmr/*solve*.sfmr"))[0]
    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    idx, dst = KdForest(desc, preset=preset).query(desc, k=K)
    idx = idx.astype(np.int64)
    img = bank.image_label.astype(np.int64)
    point = bank.point_label
    pos = point >= 0

    nb_idx = idx[:, 1:]  # 16 neighbours (drop self)
    nb_dist = dst[:, 1:]
    d1 = dst[:, 1]
    d2 = dst[:, 2]
    ratio = nb_dist / np.maximum(d2[:, None], 1e-6)  # each neighbour as ×d2

    # Only cross-image neighbours matter (clustering drops same-image ones).
    # Co-observation = same 3-D point; background = cross-image, different point.
    cross = img[nb_idx] != img[:, None]
    same_pt = (point[nb_idx] == point[:, None]) & pos[:, None]
    cotrack = cross & same_pt
    bg = cross & ~same_pt

    # In-track sources that actually have ≥1 co-observation among the 16.
    src = pos & cotrack.any(1)
    co_r = ratio[src][cotrack[src]]  # co-obs distances (×d2)
    bg_r = ratio[src][bg[src]]  # cross-image background distances (×d2)
    d2d1 = (d2 / np.maximum(d1, 1e-6))[src]

    return dict(
        name=ws.replace("_ws", ""),
        n=bank.n,
        n_src=int(src.sum()),
        d2d1=d2d1,
        co_r=co_r,
        bg_r=bg_r,
    )


def print_table(r: dict) -> None:
    co, bg = r["co_r"], r["bg_r"]
    print(f"\n{r['name']}: n={r['n']}  in-track-with-coobs={r['n_src']}")
    print(
        f"  d2/d1: median={np.median(r['d2d1']):.2f}  "
        f"p90={np.percentile(r['d2d1'], 90):.2f}"
    )
    print(
        f"  co-obs  dist/d2: p50={np.median(co):.2f}  p90={np.percentile(co, 90):.2f}"
        f"  p95={np.percentile(co, 95):.2f}"
    )
    print(
        f"  backgrnd dist/d2: p5={np.percentile(bg, 5):.2f}  "
        f"p10={np.percentile(bg, 10):.2f}  p50={np.median(bg):.2f}"
    )
    print(f"  {'c (=radius/d2)':>14} {'co-obs kept':>12} {'bg admitted':>12}")
    for c in CS:
        print(f"  {c:>14.2f} {(co <= c).mean():>11.1%} {(bg <= c).mean():>12.1%}")


def plot(results, outpath):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, r in zip(axes.ravel(), results):
        bins = np.linspace(0, 5, 60)
        ax.hist(
            r["co_r"],
            bins=bins,
            density=True,
            alpha=0.6,
            color="C2",
            label="co-observation",
        )
        ax.hist(
            r["bg_r"],
            bins=bins,
            density=True,
            alpha=0.6,
            color="C1",
            label="background",
        )
        for c in (1.5, 2.0, 3.0):
            ax.axvline(c, color="k", ls=":", alpha=0.5)
        ax.set_title(r["name"])
        ax.set_xlabel("neighbour distance / d2")
        ax.set_ylabel("density")
        ax.set_xlim(0, 5)
        ax.grid(alpha=0.3)
    axes[0, 0].legend()
    fig.suptitle(
        "Neighbour distance (÷ d2) for real tracks: co-obs vs background", y=0.997
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outpath, dpi=130)
    print(f"\nwrote {outpath}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="accurate")
    ap.add_argument("--datasets", nargs="*", default=DATASETS)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    print(
        "Where do true co-observations sit, as a multiple of d2?\n"
        "radius = c·d2 keeps a neighbour iff its distance ≤ c·d2."
    )
    results = [run_one(ws, args.preset) for ws in args.datasets]
    for r in results:
        print_table(r)
    if len(results) == len(DATASETS):
        out = Path(args.outdir)
        out.mkdir(parents=True, exist_ok=True)
        plot(results, str(out / "exp15_track_ratios.png"))


if __name__ == "__main__":
    main()
