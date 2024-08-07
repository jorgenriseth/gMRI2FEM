from pathlib import Path
from typing import Optional

import nibabel
import numpy as np
import skimage
from pint import UnitRegistry

from gmri2fem.analysis import seg_groups
from gmri2fem.utils import nan_filter_gaussian

# Define units
ureg = UnitRegistry()
ms = ureg.ms


ANSORGE_T1_TIMES = {
    "csf": 3700.0 * ms,
    "white-matter": 1084.0 * ms,
    "gray-matter": 1274.0 * ms,
    "brainstem": 993.0 * ms,  # Uncertainty: labeled spinal cord in source.
}


def label_mask(aseg: np.ndarray, labels: list[int]) -> bool:
    return sum([aseg == i for i in labels]) == 1


def brain_mask(aseg: np.ndarray, refine: bool = False) -> np.ndarray:
    mask = np.logical_and(~np.isnan(aseg), aseg != 0)
    mask = remove_csf(mask, aseg, None, replacement_value=False)
    if refine:
        mask = skimage.morphology.remove_small_objects(mask, 50, connectivity=2)
        mask = skimage.morphology.remove_small_holes(mask, 5, connectivity=2)
    return mask


def remove_csf(
    data: np.ndarray,
    aseg: np.ndarray,
    csf_labels: Optional[list[int]] = None,
    replacement_value: float | bool = np.nan,
):
    if csf_labels is None:
        csf_labels = seg_groups.default_segmentation_groups()["csf"]
    csfmask = label_mask(aseg, csf_labels)
    return np.where(csfmask, replacement_value, data)


def create_t1_map(
    aseg: np.ndarray,
    times: dict[str, float],
    default_t1: float,
    dropcsf: bool = True,
    labels: Optional[dict[str, list[int]]] = None,
) -> np.ndarray:
    if labels is None:
        labels = {
            key: val
            for key, val in seg_groups.default_segmentation_groups().items()
            if key in times
        }
    mask = brain_mask(aseg, refine=True)
    im = np.where(mask, default_t1, np.nan)
    for label, t1 in times.items():
        mask = label_mask(aseg, labels[label])
        im = np.where(mask, t1, im)

    if dropcsf:
        im = remove_csf(im, aseg, csf_labels=labels["csf"])
    return im


def replace_nan(data: np.ndarray, value: float | int | bool) -> np.ndarray:
    return np.where(np.isnan(data), value, data)


def mask_image(aseg_file: Path, output: Path) -> Path:
    aseg = nibabel.freesurfer.mghformat.load(aseg_file)
    data, affine = aseg.get_fdata(), aseg.affine
    mask = brain_mask(data, refine=True)
    mask = remove_csf(mask, data)
    mask = replace_nan(mask, False)
    np.save(output.with_suffix(".npy"), mask.astype(bool))
    nibabel.freesurfer.mghformat.save(
        nibabel.freesurfer.mghformat.MGHImage(mask.astype(np.single), affine), output
    )
    return Path(output)


def t1map_image(aseg_file: Path, output: Path, sigma: float) -> Path:
    times = {key: val.to("s").magnitude for key, val in ANSORGE_T1_TIMES.items()}
    aseg = nibabel.load(aseg_file)
    data = create_t1_map(aseg.get_fdata(), times, times["white-matter"], dropcsf=True)
    if sigma > 0:
        dilated = nan_filter_gaussian(data, sigma)
        data = np.where(
            np.logical_and(dilated > 0.0, np.isnan(data)),
            dilated,
            data,
        )
    t1map = nibabel.freesurfer.mghformat.MGHImage(
        dataobj=data.astype(np.single), affine=aseg.affine
    )
    nibabel.freesurfer.mghformat.save(t1map, output)
    return Path(output)


def create_segmentation(aseg):
    seg_numbering = {"csf": -1, "white-matter": 1, "gray-matter": 2, "brainstem": 3}
    return create_t1_map(aseg, seg_numbering, default_t1=0)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--aseg", help="Path to aseg.mgz-file.")
    parser.add_argument("--t1map", help="Path to output T1-map mgz-file")
    parser.add_argument("--mask", help="Paht to output brainmask mgz-file.")
    args = parser.parse_args()

    t1map_image(aseg_file=Path(args.aseg), output=Path(args.t1map), sigma=0)
    mask_image(aseg_file=Path(args.aseg), output=Path(args.mask))
