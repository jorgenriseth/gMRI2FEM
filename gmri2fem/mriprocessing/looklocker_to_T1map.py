import numpy as np
import scipy
import re
from pathlib import Path
from typing import Callable


import nibabel
import time
from loguru import logger


def inversion_recovery_curve(
    t: np.floating, a: float, b: float, T1: np.floating
) -> np.floating:
    return a * (1 - np.exp(-t / T1)) + b


def ll_correction(a, b, T1_pre):
    return (a / (a + b) - 1) * T1_pre


def f(t, a, b, T1):
    IR = inversion_recovery_curve(t, a, b, T1)
    y = np.abs(IR)
    return y


def voxel_fitter(t: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    def fit_voxel(TI: np.ndarray) -> np.ndarray:
        amin = max(1, min(TI.argmin(), 11))
        t_min = (0.25 * t[:-2] + 0.5 * t[1:-1] + 0.25 * t[2:])[amin - 1]
        w = 1.1
        a0 = TI[0] + w * TI[-1]
        b0 = -TI[0] * 0.9

        T10 = (-a0 / b0) * t_min
        p0 = np.array((a0, b0, T10))
        try:
            popt, pcov = scipy.optimize.curve_fit(f, xdata=t, ydata=TI, p0=p0)
            aopt, bopt, t1opt = popt
            if t1opt > t[-1] or t1opt < 0 or aopt < 0:
                return np.nan * np.ones_like(p0)
        except RuntimeError:
            return np.nan * np.ones_like(p0)
        return popt

    return fit_voxel


def estimate_t1map(ll_dir: Path, mask_quantile=0.6) -> nibabel.nifti1.Nifti1Image:
    first_p = next(ll_dir.glob("sub-*_ses-*_LookLocker_t*.nii.gz"))
    subject_id = re.match(
        r"(sub-\w+)_ses-\w+_LookLocker_t\d*.nii.gz", first_p.name
    ).groups()[0]
    session_id = re.match(
        r"sub-\w+_(ses-\w+)_LookLocker_t\d*.nii.gz", first_p.name
    ).groups()[0]
    t_labels = sorted(
        [
            float(
                re.match("sub-\w+_ses-\w+_LookLocker_t(\d*).nii.gz", p.name).groups()[0]
            )
            for p in ll_dir.glob("sub-*_ses-*_LookLocker_t*.nii.gz")
        ]
    )
    assert len(t_labels) > 0
    t_data = np.array(t_labels) / 1000.0
    lls = [
        ll_dir / f"{subject_id}_{session_id}_LookLocker_t{int(ti)}.nii.gz"
        for ti in t_labels
    ]
    mris = [nibabel.nifti1.load(ll) for ll in lls]
    data = [mri.get_fdata(dtype=np.single) for mri in mris]
    D = np.array(data)
    q = np.quantile(D[0], mask_quantile)
    valid_voxels = D[0] > q
    voxel_mask = np.where(valid_voxels)
    voxel_mask = np.array(voxel_mask).T
    Dmasked = np.array([D[:, i, j, k] for (i, j, k) in voxel_mask])
    fit_voxel = voxel_fitter(t_data)
    vfunc = np.vectorize(fit_voxel, signature="(n) -> (3)")
    tic = time.time()
    fitted_coefficients = vfunc(Dmasked)
    toc = time.time()
    logger.info(
        f"Fitted LL T1map in {ll_dir}, {len(voxel_mask)} voxel over {(toc - tic) / 60:.3f}min "
    )
    T1map = np.nan * np.zeros_like(D[0])
    T1map[voxel_mask[:, 0], voxel_mask[:, 1], voxel_mask[:, 2]] = (
        ll_correction(
            fitted_coefficients[:, 0],
            fitted_coefficients[:, 1],
            fitted_coefficients[:, 2],
        )
        * 1000.0  # convert to ms
    )
    return nibabel.nifti1.Nifti1Image(T1map, mris[0].affine)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--inputdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mask_quantile", type=float, default=0.0)
    args = parser.parse_args()
    T1map_nii = estimate_t1map(args.inputdir, args.mask_quantile)
    nibabel.nifti1.save(T1map_nii, args.output)

    if "T1map" in str(args.output):
        R1map_nii = nibabel.nifti1.Nifti1Image(
            1.0 / T1map_nii.get_fdata(dtype=np.single), T1map_nii.affine
        )
        nibabel.nifti1.save(R1map_nii, str(args.output).replace("T1map", "R1map"))
