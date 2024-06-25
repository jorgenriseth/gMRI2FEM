from pathlib import Path

import nibabel.nifti1 as nifti1
import numpy as np
import skimage


def largest_island(mask: np.ndarray, connectivity: int = 1) -> np.ndarray:
    newmask = skimage.measure.label(mask, connectivity=connectivity)
    regions = skimage.measure.regionprops(newmask)
    regions.sort(key=lambda x: x.num_pixels, reverse=True)
    return newmask == regions[0].label


def create_csf_mask(
    se0_path: Path,
    cutoff_percent: float = 25,
) -> nifti1.Nifti1Image:
    se_mri = nifti1.load(se0_path)
    se = se_mri.get_fdata()
    cutoff = (cutoff_percent / 100.0) * np.quantile(se, 0.999)
    csf_mask = se > cutoff
    csf_mask = largest_island(csf_mask)
    csf_mask = skimage.morphology.remove_small_holes(csf_mask)
    csf_mask = skimage.morphology.binary_erosion(csf_mask)
    csf_mask_nii = nifti1.Nifti1Image(csf_mask, se_mri.affine, header=se_mri.header)
    return csf_mask_nii


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--se", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    nii = create_csf_mask(args.se, cutoff_percent=25)
    nifti1.save(nii, args.output)
