# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 19 — is the "background" really background, or missed tracks?

exp17's separability ceiling and exp18's background-admission are both judged
against the reference solve's tracks. But the cluster matcher's own
(geometrically verified) solve recovers far more tracks, so the reference is
incomplete: some neighbours it labels "background" are real co-observations it
never tracked. This quantifies that.

Both solves index the *same* descriptors (same workspace .sift files); only the
per-feature track id differs. For each within-radius cross-image neighbour edge,
we compare its label under the reference solve vs under the larger solve. The key
number: of edges the reference calls background, what fraction the larger solve
calls a co-observation.

Usage:
    pixi run -e experiments python experiments/exp19_reference_incompleteness.py \
        ../dino_dog_toy_ws/sfmr/<reference>.sfmr ../dino_dog_toy_ws/sfmr/<larger>.sfmr
"""

from __future__ import annotations

import argparse

import numpy as np

from exp03_radius_clusters import K
from exp05_cluster_match import cliff_threshold, derive_threshold
from sfm_descriptors import load_descriptor_bank


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reference", help="lean reference solve .sfmr")
    ap.add_argument("larger", help="larger (cluster-matcher) solve .sfmr")
    ap.add_argument("--preset", default="accurate")
    args = ap.parse_args()
    from sfmtool import KdForest

    ref = load_descriptor_bank(args.reference)
    big = load_descriptor_bank(args.larger)
    assert ref.n == big.n, "banks must index the same descriptors"
    assert np.array_equal(ref.image_label, big.image_label)
    assert np.array_equal(ref.feature_label, big.feature_label)

    desc = np.ascontiguousarray(ref.descriptors)
    idx, dst = KdForest(desc, preset=args.preset).query(desc, k=K)
    idx = idx.astype(np.int64)
    img = ref.image_label.astype(np.int64)
    p_ref, p_big = ref.point_label, big.point_label

    print(
        f"reference tracks: {int((p_ref >= 0).sum())} in-track feats, "
        f"{len(np.unique(p_ref[p_ref >= 0]))} points"
    )
    print(
        f"larger    tracks: {int((p_big >= 0).sum())} in-track feats, "
        f"{len(np.unique(p_big[p_big >= 0]))} points"
    )

    i = np.repeat(np.arange(ref.n), K - 1)
    j = idx[:, 1:].ravel()
    d = dst[:, 1:].ravel()
    cross = img[i] != img[j]

    ref_co = cross & (p_ref[j] == p_ref[i]) & (p_ref[i] >= 0)
    big_co = cross & (p_big[j] == p_big[i]) & (p_big[i] >= 0)

    for name, T in [
        ("Otsu(d1)", derive_threshold(dst[:, 1], "otsu")),
        ("cliff p50", cliff_threshold(dst)),
    ]:
        keep = cross & (d <= T)
        kept = int(keep.sum())
        refbg = keep & ~ref_co  # within radius, reference says NOT same track
        # of those, how many does the larger solve call a co-observation
        rescued = int((refbg & big_co).sum())
        print(f"\nradius = {name} ({T:.0f}):  kept cross-image edges = {kept}")
        print(f"  reference co-obs among kept: {int((keep & ref_co).sum()) / kept:.1%}")
        print(f"  reference 'background' among kept: {int(refbg.sum()) / kept:.1%}")
        print(
            f"  ...of that 'background', larger solve calls co-obs: "
            f"{rescued / max(int(refbg.sum()), 1):.1%}  ({rescued} edges)"
        )
        union_co = int((keep & (ref_co | big_co)).sum())
        print(f"  co-obs by either solve among kept: {union_co / kept:.1%}")


if __name__ == "__main__":
    main()
