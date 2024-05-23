import re
import warnings
from functools import partial
from pathlib import Path
from typing import Callable, Any

import nibabel
import numpy as np
import scipy
import skimage
import tqdm
from loguru import logger
from scipy.optimize import OptimizeWarning

from gmri2fem.utils import nan_filter_gaussian


def f(t, x1, x2, x3):
    return np.abs(x1 * (1.0 - (1.1 + x2**2) * np.exp(-(x3**2) * t)))


@np.errstate(divide="raise", invalid="raise", over="raise")
def curve_fit_wrapper(f, t, y, p0):
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        popt, _ = scipy.optimize.curve_fit(f, xdata=t, ydata=y, p0=p0, maxfev=1000)
    return popt


def fit_voxel(time_s: np.ndarray, pbar, m: np.ndarray) -> np.ndarray:
    if pbar is not None:
        pbar.update(1)
    p0 = np.array((1.0, 1.0, 1.0))
    if not np.all(np.isfinite(m)):
        return np.nan * np.zeros_like(p0)
    try:
        popt = curve_fit_wrapper(f, time_s, m, p0)
    except (OptimizeWarning, FloatingPointError):
        return np.nan * np.zeros_like(p0)
    except RuntimeError as e:
        if "maxfev" in str(e):
            return np.nan * np.zeros_like(p0)
        raise e
    return popt


def parse_looklocker_path(
    p: Path, raise_error: bool = True
) -> tuple[str | Any, ...] | None:
    pattern = r"sub-(\w+)_ses-(\w+)_LookLocker_t(\d+).nii.gz"
    m = re.match(pattern, str(p.name))
    if m is not None:
        return m.groups()
    if raise_error:
        raise ValueError(f"Pattern {pattern} not found in {p}")
    return None


def estimate_t1map(
    ll_dir: Path, threshold_number: int = 0
) -> nibabel.nifti1.Nifti1Image:
    match_groups = [
        mgroup
        for mgroup in map(
            parse_looklocker_path, ll_dir.glob("sub-*_ses-*_LookLocker_t*.nii.gz")
        )
        if mgroup is not None
    ]
    sub, ses = match_groups[0][0], match_groups[0][1]
    t_labels = sorted([float(g[2]) for g in match_groups])
    t_data = np.array(t_labels) / 1000.0
    lls = [
        ll_dir / f"sub-{sub}_ses-{ses}_LookLocker_t{int(ti)}.nii.gz" for ti in t_labels
    ]
    D = np.array([nibabel.nifti1.load(ll).get_fdata() for ll in lls])
    lowermost = np.sort(np.unique(D[0]))

    # Keep only largest island, but with hole-filling.
    mask_labels = skimage.measure.label(D[0] > lowermost[threshold_number])
    regions = skimage.measure.regionprops(mask_labels)
    regions.sort(key=lambda x: x.num_pixels, reverse=True)
    largest_island = mask_labels == regions[0].label
    largest_island = skimage.morphology.remove_small_holes(
        largest_island, 10 ** (D.ndim - 1), connectivity=2
    )
    valid_voxels = (D[0] > 0) * largest_island

    D_normalized = np.nan * np.zeros_like(D)
    D_normalized[:, valid_voxels] = D[:, valid_voxels] / D[0][valid_voxels]
    voxel_mask = np.array(np.where(valid_voxels)).T
    Dmasked = np.array([D_normalized[:, i, j, k] for (i, j, k) in voxel_mask])

    with tqdm.tqdm(total=len(Dmasked)) as pbar:
        voxel_fitter = partial(fit_voxel, t_data, pbar)
        vfunc = np.vectorize(voxel_fitter, signature="(n) -> (3)")
        fitted_coefficients = vfunc(Dmasked)

    x1, x2, x3 = (
        fitted_coefficients[:, 0],
        fitted_coefficients[:, 1],
        fitted_coefficients[:, 2],
    )

    I, J, K = voxel_mask.T
    T1map = np.nan * np.zeros_like(D[0])
    T1map[I, J, K] = (0.1 + x2**2) / x3**2 * 1000.0  # convert to ms

    affine = nibabel.nifti1.load(lls[0]).affine
    return nibabel.nifti1.Nifti1Image(T1map.astype(np.single), affine)


def postprocess_T1map(
    T1map_mri: nibabel.nifti1.Nifti1Image,
    T1_lo: float,
    T1_hi: float,
    radius: int = 10,
    erode_dilate_factor: float = 1.3,
) -> nibabel.nifti1.Nifti1Image:
    T1map = T1map_mri.get_fdata()

    # Create mask for largest island.
    mask = skimage.measure.label(np.isfinite(T1map))
    regions = skimage.measure.regionprops(mask)
    regions.sort(key=lambda x: x.num_pixels, reverse=True)
    mask = mask == regions[0].label
    skimage.morphology.remove_small_holes(
        mask, 10 ** (mask.ndim), connectivity=2, out=mask
    )
    skimage.morphology.binary_dilation(mask, skimage.morphology.ball(radius), out=mask)
    skimage.morphology.binary_erosion(
        mask, skimage.morphology.ball(erode_dilate_factor * radius), out=mask
    )
    # Remove outliers and small surface artifacts
    surface_vox = np.isfinite(T1map) * (~mask)
    print(f"Removing {surface_vox.sum()} voxels outside of the head mask")
    T1map[~mask] = np.nan

    outliers = np.logical_or(T1map < T1_lo, T1_hi < T1map)
    print("Removing", outliers.sum(), f"voxels outside the range ({T1_lo}, {T1_hi}).")
    T1map[outliers] = np.nan

    # Fill internallly missing values
    fill_mask = np.isnan(T1map) * mask
    print(f"Filling in {fill_mask.sum()} voxels within the mask.")
    T1map[fill_mask] = nan_filter_gaussian(T1map, 1.0)[fill_mask]
    return nibabel.nifti1.Nifti1Image(T1map, T1map_mri.affine)


def T1_to_R1(T1map_mri: nibabel.nifti1.Nifti1Image):
    T1map = T1map_mri.get_fdata()
    return nibabel.nifti1.Nifti1Image(1 / T1map, T1map_mri.affine)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--inputdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threshold_number", type=int, default=0)
    parser.add_argument("--T1_low", type=int, default=50)
    parser.add_argument("--T1_hi", type=int, default=5000)
    args = parser.parse_args()

    T1map_nii = estimate_t1map(args.inputdir, args.threshold_number)
    nibabel.nifti1.save(T1map_nii, args.output)

    postprocess_T1map(
        T1map_nii, args.T1_low, args.T1_hi, radius=10, erode_dilate_factor=1.3
    )
