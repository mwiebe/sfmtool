# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 04 — derive the clustering radius from the data (no labels).

exp03 hand-picked T≈80. But the nearest-neighbour distance d1 is bimodal —
a "has-a-real-neighbour" mode and an "isolated" mode — so the antimode between
them is a label-free threshold. We pick it two ways (Otsu on the d1 histogram,
and the crossover of a 2-component Gaussian mixture on log d1), then validate the
resulting clustering against the solve's tracks.

Uses the 17-NN cache written by exp03.

Usage:
    pixi run -e experiments python experiments/exp04_auto_threshold.py \
        ../dino_dog_toy_ws/sfmr/*.sfmr
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.mixture import GaussianMixture

from exp03_radius_clusters import K, build_and_score, knn_all
from sfm_descriptors import load_descriptor_bank


def otsu_threshold(x: np.ndarray, nbins: int = 256) -> float:
    """Otsu's between-class-variance threshold on a 1-D sample."""
    hist, edges = np.histogram(x, bins=nbins)
    centers = (edges[:-1] + edges[1:]) / 2
    w = hist.astype(float) / hist.sum()
    cumw = np.cumsum(w)
    cumm = np.cumsum(w * centers)
    gm = cumm[-1]
    denom = cumw * (1 - cumw)
    num = (gm * cumw - cumm) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_b = np.where(denom > 0, num / denom, 0.0)
    return float(centers[np.argmax(sigma_b)])


def gmm_crossover(x: np.ndarray) -> float:
    """Threshold where a 2-component GMM (on log x) switches dominant component."""
    lx = np.log(np.clip(x, 1e-3, None)).reshape(-1, 1)
    gm = GaussianMixture(n_components=2, random_state=0).fit(lx)
    lo, hi = sorted(gm.means_.ravel())
    grid = np.linspace(lo, hi, 2000).reshape(-1, 1)
    pred = gm.predict(grid).ravel()
    # first grid point assigned to the high-mean component
    hi_comp = np.argmax(gm.means_.ravel())
    cross = grid[pred == hi_comp]
    return float(np.exp(cross[0, 0])) if len(cross) else float(np.exp((lo + hi) / 2))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("sfmr")
    ap.add_argument("--cache", default=None)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    path = sorted(glob.glob(args.sfmr))[0]
    bank = load_descriptor_bank(path)
    cache = Path(args.cache or f"out/{Path(bank.workspace_dir).name}_knn17.npz")
    if cache.exists():
        z = np.load(cache)
        idx, dst = z["idx"], z["dst"]
    else:
        print(f"computing exact {K}-NN over {bank.n} (one-time) -> {cache}")
        idx, dst = knn_all(bank.descriptors, K)
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez(cache, idx=idx, dst=dst)
    print(f"{Path(bank.workspace_dir).name}: n={bank.n} "
          f"in_track={int(bank.in_track.sum())}")

    d1 = dst[:, 1]
    d5 = dst[:, 5]
    # Fit on d1 excluding exact-duplicate spike (d1~0) and the long isolated tail,
    # so neither artefact captures a mixture component.
    d1f = d1[(d1 > 1.0) & (d1 <= np.percentile(d1, 99.5))]

    t_otsu = otsu_threshold(d1f)
    t_gmm = gmm_crossover(d1f)
    print(f"d1 distribution: median={np.median(d1):.1f} "
          f"p90={np.percentile(d1, 90):.1f}")
    print(f"derived thresholds (label-free):  Otsu={t_otsu:.1f}  GMM={t_gmm:.1f}")

    # data-derived prefilter too: isolated = d1 above the antimode, or flat ratio.
    ratio = d1 / np.maximum(d5, 1e-6)
    t_pf = t_gmm  # reuse the antimode as the isolated cutoff
    active = ~((ratio > 0.85) | (d1 > t_pf))
    print(f"prefilter (d1>{t_pf:.0f} or ratio>0.85) keeps "
          f"{active.mean():.1%}")

    print("\nclustering at the derived thresholds (with prefilter):")
    print(f"  {'T(source)':>12} {'comps':>7} {'purity':>7} {'recov':>7} "
          f"{'imgUniq':>8} {'cover':>7}")
    for name, T in [("Otsu", t_otsu), ("GMM", t_gmm), ("hand=80", 80.0)]:
        r = build_and_score(bank, idx, dst, T, active)
        print(f"  {name + f'={T:.0f}':>12} {r.get('components', 0):>7} "
              f"{r.get('dom_cov', 0):>7.3f} {r.get('recovery', 0):>7.3f} "
              f"{r.get('img_unique_frac', 0):>8.1%} "
              f"{r.get('in_track_coverage', 0):>7.1%}")

    # plot d1 with derived thresholds, split by ground truth for validation
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    in_track = bank.point_label >= 0
    plt.figure(figsize=(8, 5))
    bins = np.linspace(0, np.percentile(d1, 99), 80)
    plt.hist(d1[in_track], bins=bins, density=True, alpha=0.55, label="in-track")
    plt.hist(d1[~in_track], bins=bins, density=True, alpha=0.55, label="background")
    plt.axvline(t_otsu, color="k", ls="--", label=f"Otsu={t_otsu:.0f}")
    plt.axvline(t_gmm, color="r", ls=":", label=f"GMM={t_gmm:.0f}")
    plt.xlabel("nearest-neighbour distance d1 (L2)")
    plt.ylabel("density")
    name = Path(bank.workspace_dir).name
    plt.title(f"{name}: data-derived threshold from the bimodal d1 distribution")
    plt.legend()
    plt.tight_layout()
    f = outdir / f"exp04_{name}_threshold.png"
    plt.savefig(f, dpi=110)
    plt.close()
    print(f"\nwrote {f}")


if __name__ == "__main__":
    main()
