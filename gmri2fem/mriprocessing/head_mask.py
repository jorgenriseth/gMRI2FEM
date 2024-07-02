import numpy as np
import nibabel.nifti1 as nifti1

from pathlib import Path
from gmri2fem.utils import mri_facemask


if __name__ == "__main__"
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

# args.input = Path("data/mri_dataset/sub-01/ses-01/anat/sub-01_ses-01_T1w.nii.gz")
# args.output = Path("testmask.nii.gz")

    nii = nifti1.load(args.input)
    vol = nii.get_fdata(dtype=np.single)
    binary = mri_facemask(vol)

    mask_nii = nifti1.Nifti1Image(binary.astype(np.uint8), nii.affine)
    nifti1.save(mask_nii, args.output)
