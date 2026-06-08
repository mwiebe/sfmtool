# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Shared loader for descriptor-space matching experiments.

Pulls every SIFT descriptor from every image in a workspace into one big
``(N, 128)`` array, and labels each descriptor with the solve's ground truth:
the 3D point id of the track it belongs to (or -1 if the solve never used it).

This is deliberately standalone POC plumbing — it lives under ``experiments/``
and is not wired into the ``sfmtool`` package.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from sfmtool._sfmtool import SfmrReconstruction
from sfmtool.sift.file import SiftReader, get_sift_path_for_image


@dataclass
class DescriptorBank:
    """All descriptors from a reconstruction, flattened and labelled.

    Attributes:
        descriptors: ``(N, 128)`` uint8 array, every feature from every image,
            concatenated image-by-image.
        image_label: ``(N,)`` int32, the image index each descriptor came from.
        feature_label: ``(N,)`` int32, the per-image feature index (row in that
            image's ``.sift`` file).
        point_label: ``(N,)`` int64, the solve's 3D point id for the track this
            descriptor belongs to, or -1 if the descriptor was never used.
        image_names: list of image names (relative to the workspace).
        workspace_dir: workspace root path.
    """

    descriptors: np.ndarray
    image_label: np.ndarray
    feature_label: np.ndarray
    point_label: np.ndarray
    image_names: list[str]
    workspace_dir: Path

    @property
    def n(self) -> int:
        return len(self.descriptors)

    @property
    def in_track(self) -> np.ndarray:
        """Boolean mask of descriptors the solve actually used in a track."""
        return self.point_label >= 0


def load_descriptor_bank(sfmr_path: str | Path) -> DescriptorBank:
    """Load every descriptor in the reconstruction's images, labelled by track.

    Args:
        sfmr_path: path to the ``.sfmr`` reconstruction.

    Returns:
        A populated :class:`DescriptorBank`.
    """
    recon = SfmrReconstruction.load(str(sfmr_path))
    ws = Path(recon.workspace_dir)
    names = list(recon.image_names)

    # Solve ground truth: parallel observation arrays.
    track_img = np.asarray(recon.track_image_indexes, dtype=np.int64)
    track_feat = np.asarray(recon.track_feature_indexes, dtype=np.int64)
    track_pid = np.asarray(recon.track_point_ids, dtype=np.int64)

    desc_blocks: list[np.ndarray] = []
    image_label_blocks: list[np.ndarray] = []
    feature_label_blocks: list[np.ndarray] = []
    # (image_idx, feat_idx) -> global row, so we can scatter track labels.
    base_offset: list[int] = []

    offset = 0
    for img_idx, name in enumerate(names):
        sift_path = get_sift_path_for_image(str(ws / name))
        with SiftReader(sift_path) as reader:
            desc = np.asarray(reader.read_descriptors())
        k = len(desc)
        desc_blocks.append(desc)
        image_label_blocks.append(np.full(k, img_idx, dtype=np.int32))
        feature_label_blocks.append(np.arange(k, dtype=np.int32))
        base_offset.append(offset)
        offset += k

    descriptors = np.concatenate(desc_blocks, axis=0)
    image_label = np.concatenate(image_label_blocks, axis=0)
    feature_label = np.concatenate(feature_label_blocks, axis=0)
    point_label = np.full(offset, -1, dtype=np.int64)

    base_offset_arr = np.asarray(base_offset, dtype=np.int64)
    global_rows = base_offset_arr[track_img] + track_feat
    point_label[global_rows] = track_pid

    return DescriptorBank(
        descriptors=descriptors,
        image_label=image_label,
        feature_label=feature_label,
        point_label=point_label,
        image_names=names,
        workspace_dir=ws,
    )
