# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Compare the per-point floor radius to the real track extent in descriptor space.

For every in-track descriptor we measure two distances:
  - floor radius  = alpha * dist[d]  (the membership radius the rule uses),
  - track extent  = distance to its nearest / farthest same-track co-observation
    (how far the radius must reach to capture one / all co-observations).

Overlays their distributions so we can see whether the floor is sized like the
tracks or runs generous. Writes a small PNG.

Usage:
    pixi run -e experiments python experiments/make_radius_plot.py \
        --solve ../seattle_backyard_ws/sfmr/*solve*.sfmr --out /tmp/radii.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from sfm_descriptors import load_descriptor_bank, resolve_solve
from sfmtool import KdForest

D_RANK = 28
BG_ALPHA = 0.8


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--solve", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    bank = load_descriptor_bank(resolve_solve(args.solve))
    desc = np.ascontiguousarray(bank.descriptors)
    _, dst = KdForest(desc, preset="accurate").query(desc, k=D_RANK + 4)
    floor = BG_ALPHA * dst[:, D_RANK]

    pt = bank.point_label
    Xf = desc.astype(np.float32)
    near = np.full(bank.n, np.nan, np.float32)
    far = np.full(bank.n, np.nan, np.float32)
    order = np.argsort(pt, kind="stable")
    pts = pt[order]
    bounds = np.r_[0, np.flatnonzero(np.diff(pts)) + 1, len(pts)]
    for s, e in zip(bounds[:-1], bounds[1:]):
        if pts[s] < 0 or e - s < 2:
            continue
        members = order[s:e]
        D = np.linalg.norm(Xf[members][:, None, :] - Xf[members][None, :, :], axis=2)
        np.fill_diagonal(D, np.inf)
        near[members] = D.min(1)
        D[~np.isfinite(D)] = -np.inf
        far[members] = D.max(1)

    m = np.isfinite(near)
    bins = np.linspace(0, np.percentile(np.r_[floor[m], far[m]], 99), 70)
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    for vals, color, lbl in [
        (near[m], "#2a9d8f", "nearest co-observation"),
        (far[m], "#e76f51", "farthest co-observation"),
        (floor[m], "#264653", "floor radius α·B"),
    ]:
        ax.hist(np.clip(vals, 0, bins[-1]), bins=bins, histtype="step", lw=1.8,
                color=color, label=f"{lbl} (median {np.median(vals):.0f})")
    ax.set_xlabel("descriptor distance")
    ax.set_ylabel("in-track descriptors")
    ax.set_xlim(0, bins[-1])
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(
        f"wrote {out}  floor_med={np.median(floor[m]):.0f} "
        f"near_med={np.median(near[m]):.0f} far_med={np.median(far[m]):.0f} "
        f"captured_near={np.mean(floor[m] >= near[m]):.1%} "
        f"captured_far={np.mean(floor[m] >= far[m]):.1%}"
    )


if __name__ == "__main__":
    main()
