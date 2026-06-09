# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Approach-section figure: the background floor between co-observations and
background, on a few example descriptors from one reconstruction.

For descriptors drawn from tracks of a few sizes, plot the sorted neighbour
distances by rank, colouring co-observations vs background, and mark the
background scale B (the d-th nearest) and the membership radius alpha*B. Writes a
small palette-quantised PNG for specs/core/images/.

Usage:
    pixi run -e experiments python experiments/make_floor_profile.py \
        --solve ../seattle_backyard_ws/sfmr/*solve*.sfmr \
        --out ../specs/core/images/floor-profile.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from PIL import Image

from sfm_descriptors import load_descriptor_bank, resolve_solve
from sfmtool import KdForest

D_RANK = 28
BG_ALPHA = 0.8
KQ = 44
GREEN, RED, DARK, GREY = "#2a9d8f", "#e76f51", "#264653", "#b9bec6"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--solve", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sizes", type=int, nargs="*", default=[3, 8, 16])
    args = ap.parse_args()

    bank = load_descriptor_bank(resolve_solve(args.solve))
    desc = np.ascontiguousarray(bank.descriptors)
    idx, dst = KdForest(desc, preset="accurate").query(desc, k=KQ + 1)
    idx = idx.astype(np.int64)
    pt = bank.point_label
    tsize = np.bincount(pt[pt >= 0])
    avail = np.unique(tsize[tsize >= 2])

    fig, axes = plt.subplots(1, len(args.sizes), figsize=(7.4, 2.7), sharex=True)
    ranks = np.arange(1, KQ + 1)
    for k, (ax, target) in enumerate(zip(axes, args.sizes)):
        ts = int(avail[np.argmin(np.abs(avail - target))])
        cand = np.flatnonzero((pt >= 0) & (tsize[np.clip(pt, 0, None)] == ts))
        i = int(cand[len(cand) // 2])
        dd, jj = dst[i, 1:], idx[i, 1:]
        is_co = (pt[jj] >= 0) & (pt[jj] == pt[i])
        B = dst[i, D_RANK]
        ax.scatter(ranks[~is_co], dd[~is_co], s=9, color=GREY)
        ax.scatter(ranks[is_co], dd[is_co], s=16, color=GREEN, zorder=3)
        ax.axhline(B, color=DARK, ls=":", lw=1.1)
        ax.axhline(BG_ALPHA * B, color=RED, ls="--", lw=1.3)
        ax.set_title(f"track size {ts}", fontsize=9)
        ax.set_xlabel("neighbour rank", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
        if k == 0:
            ax.set_ylabel("descriptor distance", fontsize=8)
    axes[-1].annotate(
        "background scale B", (KQ, B), fontsize=7, color=DARK, ha="right", va="bottom"
    )
    axes[-1].annotate(
        "floor α·B", (KQ, BG_ALPHA * B), fontsize=7, color=RED, ha="right", va="bottom"
    )

    handles = [
        Line2D([], [], marker="o", ls="", color=GREEN, label="co-observation"),
        Line2D([], [], marker="o", ls="", color=GREY, label="background"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=2,
        frameon=False,
        fontsize=8,
        bbox_to_anchor=(0.5, -0.02),
    )
    fig.tight_layout(rect=(0, 0.06, 1, 1))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp.png")
    fig.savefig(tmp, dpi=130)
    plt.close(fig)
    Image.open(tmp).convert("RGB").quantize(colors=64, method=Image.MAXCOVERAGE).save(
        out, optimize=True
    )
    tmp.unlink()
    print(f"wrote {out} ({out.stat().st_size // 1024} KiB)")


if __name__ == "__main__":
    main()
