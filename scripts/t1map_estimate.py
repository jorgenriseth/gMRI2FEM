import re
from pathlib import Path

import matplotlib.pyplot as plt
import nibabel
import numpy as np
import scipy


def inversion_recovery_curve(t, a, b, T1):
    return a * (1 - np.exp(-t / T1)) + b


def ll_correction(a, b, T1_pre):
    return (a / (a + b) - 1) * T1_pre


def f(t, a, b, T1):
    IR = inversion_recovery_curve(t, a, b, T1)
    return np.abs(IR)


def process_voxel(TI: np.ndarray) -> np.ndarray:
    if TI.max() < 1e-8:
        return np.nan * np.zeros(3)
    TI /= TI.max()
    t_min = t_data[TI.argmin()]
    a0 = TI[0] + TI[-1]
    b0 = -TI[0]
    T10 = -t_min / np.log(1 + b0 / a0)
    p0 = np.array((a0, b0, T10))
    try:
        popt, pcov = scipy.optimize.curve_fit(f, xdata=t_data, ydata=TI, p0=p0)
    except RuntimeError:
        return np.nan * np.zeros_like(p0)
    return popt


def bounding_box_mask(bbox: list[tuple[int, int]]) -> np.ndarray:
    bbox_voxels = np.ones(shape, dtype=bool)
    all_indices = np.indices(shape)
    for ax, lims in enumerate(bbox):
        m = (lims[0] < all_indices[ax]) * (all_indices[ax] < lims[1])
        bbox_voxels *= m
    return bbox_voxels


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    for p in args.input:
        print(p)
    t_labels = sorted(
        [float(re.match(r"(.*)_t(\d*).nii", p.name).groups()[1]) for p in args.input]
    )
    t_data = np.array(t_labels) / 1000

    raise NotImplementedError("Need to add sorting of input files")

    lls = 0
    mris = [nibabel.nifti1.load(ll) for ll in args.input]
    D = np.array([mri.get_fdata(dtype=np.single) for mri in mris])
    shape = D[0].shape

    bbox = [(10, 170), (25, 200), (70, 220)]
    bbox_voxels = bounding_box_mask(bbox)
    invalid_voxels = np.all(D == 0, axis=0) + np.any(np.isnan(D), axis=0) + ~bbox_voxels
    voxel_mask = np.where(~invalid_voxels)
    voxel_mask = np.array(voxel_mask).T
    ordered_data = np.array(D[:, i, j, k] for i, j, k in voxel_mask)
    vectorized_process_voxel = np.vectorize(process_voxel, signature="(n) -> (3)")
    print("Creating T1map. This can take a while (probably under an hour.)")
    coefficients = vectorized_process_voxel(ordered_data)

    T1map = np.nan * np.zeros_like(D[0])
    T1map[voxel_mask[:, 0], voxel_mask[:, 1], voxel_mask[:, 2]] = ll_correction(
        coefficients[:, 0], coefficients[:, 1], coefficients[:, 2]
    )

    T1map_mri = nibabel.nifti1.Nifti1Image(T1map, mris[0].affine, mris[0].header)
    nibabel.nifti1.save(T1map_mri, args.output)
