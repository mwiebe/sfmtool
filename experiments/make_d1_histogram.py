# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Render the d1 (nearest-other-descriptor distance) histogram for the spec.

Loads every descriptor of one dataset, computes `d1` for each via the KdForest,
and plots its histogram split by whether the solve used the descriptor in a track
— the "has-a-real-neighbour" mode vs the "isolated" mode. Writes a small,
palette-quantised PNG for `specs/core/images/`.

Usage:
    pixi run -e experiments python experiments/make_d1_histogram.py \
        --solve ../seattle_backyard_ws/sfmr/*solve*.sfmr \
        --out ../specs/core/images/d1-histogram-seattle.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from exp04_auto_threshold import otsu_threshold
from sfm_descriptors import load_descriptor_bank, resolve_solve
from sfmtool import KdForest


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--solve", required=True, help="path/glob to a *solve*.sfmr")
    ap.add_argument("--out", required=True, help="output .png path")
    ap.add_argument("--clip-pct", type=float, default=99.0)
    args = ap.parse_args()

    bank = load_descriptor_bank(resolve_solve(args.solve))
    desc = np.ascontiguousarray(bank.descriptors)
    _, dst = KdForest(desc, preset="accurate").query(desc, k=2)
    d1 = dst[:, 1]

    hi = float(np.percentile(d1, args.clip_pct))
    bins = np.linspace(0, hi, 80)
    t = otsu_threshold(d1[(d1 > 1.0) & (d1 <= hi)])

    fig, ax = plt.subplots(figsize=(6.0, 3.0))
    counts, edges, _ = ax.hist(
        d1, bins=bins, range=(0, hi), color="#2a6f97", edgecolor="none"
    )
    ax.axvline(t, color="#9b2226", lw=1.2, ls="--")
    ax.set_xlabel("d1 — distance to nearest other descriptor")
    ax.set_ylabel("descriptors")
    ax.set_xlim(0, hi)
    ax.set_ylim(0, counts.max() * 1.18)
    ax.spines[["top", "right"]].set_visible(False)

    # Label the two modes on their peaks, and the valley.
    mids = (edges[:-1] + edges[1:]) / 2
    left = slice(0, int(np.searchsorted(mids, t)))
    right = slice(int(np.searchsorted(mids, t)), len(mids))
    lpk = left.start + int(np.argmax(counts[left]))
    rpk = right.start + int(np.argmax(counts[right]))
    ax.annotate(
        "has a near\nneighbour",
        (mids[lpk], counts[lpk]),
        ha="center",
        va="bottom",
        fontsize=9,
        xytext=(0, 6),
        textcoords="offset points",
    )
    ax.annotate(
        "isolated",
        (mids[rpk], counts[rpk]),
        ha="center",
        va="bottom",
        fontsize=9,
        xytext=(0, 6),
        textcoords="offset points",
    )
    ax.text(
        t,
        counts.max() * 1.12,
        f"antimode ≈ {t:.0f} ",
        color="#9b2226",
        fontsize=8,
        ha="right",
        va="top",
    )
    fig.tight_layout()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp.png")
    fig.savefig(tmp, dpi=120)
    plt.close(fig)

    # Palette-quantise for a small line-art PNG.
    img = Image.open(tmp).convert("RGB").quantize(colors=32, method=Image.MAXCOVERAGE)
    img.save(out, optimize=True)
    tmp.unlink()
    print(f"wrote {out} ({out.stat().st_size // 1024} KiB), Otsu≈{t:.0f}, n={bank.n}")


if __name__ == "__main__":
    main()
