import itertools
from pathlib import Path
from typing import Literal, Sequence

import click
import numpy as np
import scipy
import skimage
import simple_mri as sm

from gmri2fem.utils import largest_island

# LUT-values of regions to look for to find slice in R,A,S-directions respectively
AXES_SEG_LABELS = {
    "left": [
        1012,
        1027,
        18,
    ],
    "right": [
        2012,
        2027,
        54,
    ],
}
ORBIT_DIMS_IN_MM = [10, 15, 10]  # along R, A, S


@click.command()
@click.option("--T1w_dir", type=Path, required=True)
@click.option("--segmentation", type=Path, required=True)
@click.option("--output", type=Path, required=True)
@click.option("--side", type=str, default="left")
def orbital_refroi(t1w_dir, segmentation, output, side):
    # FIXME: Figure out a better way to only include the
    # desired files, as the lumbar/thoracic images might
    # interfere here.
    image_paths = sorted(t1w_dir.glob("*T1w*"))
    ref_path = image_paths[0]
    ref_mri = sm.load_mri(ref_path, dtype=np.single)
    ref_affine = ref_mri.affine

    mris = np.array([sm.load_mri(path, dtype=np.single) for path in image_paths])
    seg_mri = sm.load_mri(segmentation, dtype=np.int16)

    weights = orbital_weights_distribution(ref_mri, seg_mri, side)
    masks = [threshold_weighted_image(mri.data, weights) for mri in mris]
    roi = aggregate_session_masks(masks)
    refroi_mri = sm.SimpleMRI(roi, ref_affine)
    sm.save_mri(refroi_mri, output, dtype=np.uint8)


def orbital_weights_distribution(ref_mri, seg_mri, side: Literal["left", "right"]):
    axes_seg_indices = AXES_SEG_LABELS[side]
    centers_seg = find_axial_label_centers(seg_mri, axes_seg_indices)
    centers_ras = sm.apply_affine(seg_mri.affine, centers_seg)

    ref_affine = ref_mri.affine
    std = np.array(ORBIT_DIMS_IN_MM) / 2

    # We only need the distribution in a region surrounding the
    # orbital fat, and therefore create an ealuation window to
    # avoid evaluation in the rest of the image.
    lower_bounds = sm.apply_affine(np.linalg.inv(ref_affine), centers_ras - 5 * std)
    upper_bounds = sm.apply_affine(np.linalg.inv(ref_affine), centers_ras + 5 * std)
    eval_indices = np.fromiter(
        itertools.product(
            *[
                np.arange(int(ai), int(bi))
                for ai, bi in zip(np.rint(lower_bounds), np.rint(upper_bounds))
            ]
        ),
        dtype=np.dtype((int, 3)),
    )
    G = scipy.stats.multivariate_normal(mean=centers_ras, cov=std**2)
    weights = np.zeros_like(ref_mri.data)
    I, J, K = eval_indices.T  # noqa: E741
    weights[I, J, K] = G.pdf(sm.apply_affine(ref_affine, eval_indices))
    return weights


def find_axial_label_centers(seg_mri: sm.SimpleMRI, axes_seg_labels: Sequence[int]):
    # Create weighting distribution
    centers_seg = np.array(
        [
            [x.mean() for x in np.where(seg_mri.data == label)][idx]
            for idx, label in enumerate(axes_seg_labels)
        ]
    )
    return centers_seg


def threshold_weighted_image(volume: np.ndarray, weighting: np.ndarray):
    image = volume * weighting
    (hist, bin_edges) = np.histogram(image[(image > 1e-8)], bins=256)
    bins = (bin_edges[:-1] + bin_edges[1:]) / 2
    yen_thresh = skimage.filters.threshold_yen(hist=(hist, bins))
    mask = image > yen_thresh
    if mask.sum() < 1000:
        otsu_thresh = skimage.filters.threshold_otsu(image)
        return image > otsu_thresh
    return mask


def aggregate_session_masks(masks: list[np.ndarray]) -> np.ndarray:
    roi = np.prod(masks, axis=0).astype(bool)
    roi = skimage.morphology.binary_erosion(roi, footprint=skimage.morphology.ball(1))
    # temporal_std = (
    #     vols[:, roi] / np.median(vols[:, roi], axis=1, keepdims=True)
    # ).std(axis=0)
    # roi[roi] = (
    #     np.abs(vols[0, roi] / np.median(vols[0, roi]) - 1) < 0.25
    # ) * (temporal_std < 0.1)
    return largest_island(roi)


if __name__ == "__main__":
    orbital_refroi()
