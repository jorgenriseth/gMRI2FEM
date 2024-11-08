from pathlib import Path

import nibabel.nifti1 as nifti1
import numpy as np
import skimage


def create_csf_mask(vol: np.ndarray, connectivity: int = 2, use_li: bool = False):
    if use_li:
        thresh = skimage.filters.threshold_li(vol)
        binary = vol > thresh
        binary = largest_island(binary, connectivity=2)
    else:
        (hist, bins) = np.histogram(vol[vol > 0], bins=512)
        thresh = skimage.filters.threshold_yen(hist=(hist, bins))
        binary = vol > thresh
        binary = largest_island(binary, connectivity=2)
    return binary


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--connectivity", type=int, default=2)
    parser.add_argument("--li", action="store_true")
    args = parser.parse_args()

    nii = nifti1.load(args.input)
    vol = nii.get_fdata(dtype=np.single)
    mask = create_csf_mask(vol, args.connectivity, args.li)
    mask_nii = nifti1.Nifti1Image(mask.astype(np.uint8), nii.affine)
    nifti1.save(mask_nii, args.output)
