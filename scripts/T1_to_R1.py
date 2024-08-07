from pathlib import Path
import numpy as np
import nibabel

from gmri2fem.mriprocessing.looklocker_to_T1map import T1_to_R1


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scale", type=float, default=1000)
    args = parser.parse_args()

    T1map_nii = nibabel.nifti1.load(args.input)
    R1map_nii = T1_to_R1(T1map_nii, args.scale)
    nibabel.nifti1.save(R1map_nii, args.output)
