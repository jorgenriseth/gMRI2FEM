import argparse
from pathlib import Path

import nibabel
import numpy as np


def concentration_from_T1(T1: np.ndarray, T1_0: np.ndarray, r1: float) -> np.ndarray:
    C = 1 / r1 * (1 / T1 - 1 / T1_0)
    return C


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--r1", type=float, default=3.2)
    parser.add_argument("--mask", type=Path, default=None)
    args = parser.parse_args()

    reference_volume = nibabel.freesurfer.mghformat.load(args.reference)
    reference_affine = reference_volume.header.get_vox2ras()
    reference_data = reference_volume.get_fdata(dtype=np.float32)

    volume = nibabel.freesurfer.mghformat.load(args.input)
    assert np.allclose(
        reference_affine, volume.header.get_vox2ras()
    ), "Affine transformations differ, are you sure the images are registered properly?"
    data = volume.get_fdata(dtype=np.single)

    if args.mask is not None:
        mask = nibabel.freesurfer.mghformat.load(args.mask)
        assert np.allclose(
            reference_affine, mask.header.get_vox2ras()
        ), "Affine transformations differ, are you sure the baseline and T1 Map are registered properly?"
        mask = mask.get_fdata(dtype=np.single).astype(bool)
        reference_data *= mask
        data *= mask
    else:
        mask = (reference_data > 1e-10) * (data > 1e-10)
        reference_data[~mask] = np.nan
        data[~mask] = np.nan

    concentration = concentration_from_T1(T1=data, T1_0=reference_data, r1=args.r1)
    nibabel.freesurfer.mghformat.save(
        nibabel.freesurfer.mghformat.MGHImage(concentration, reference_affine),
        args.output,
    )
