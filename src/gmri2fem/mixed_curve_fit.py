import warnings
from pathlib import Path

import nibabel
import numpy as np
import pandas as pd
import scipy
import tqdm
from loguru import logger
from scipy.optimize import OptimizeWarning


def f(t, a, b, c):
    return a * t ** (b) * np.exp(-t / c)


@np.errstate(divide="raise", invalid="raise", over="raise")
def curve_fit_wrapper(f, t, y, p0, bounds):
    """Raises error instead of catching numpy warnings, such that
    these cases may be treated."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        popt, _ = scipy.optimize.curve_fit(
            f, xdata=t, ydata=y, p0=p0, maxfev=2000, bounds=bounds
        )
    return popt, _


def fit_voxel(
    time_s: np.ndarray, pbar, m: np.ndarray, f, p0, bounds=None
) -> np.ndarray:
    if pbar is not None:
        pbar.update(1)
    if not np.all(np.isfinite(m)):
        return np.nan * np.zeros_like(p0)
    try:
        popt, _ = curve_fit_wrapper(f, time_s, m, p0, bounds=bounds)
    except (OptimizeWarning, FloatingPointError):
        return np.nan * np.zeros_like(p0)
    except RuntimeError as e:
        if "maximum number of function evaluations" in str(e):
            logger.warning("Max func-eval, setting all nan.")
            return np.nan * np.zeros_like(p0)
        raise e
    except ValueError as e:
        if "infeasible" in str(e):
            popt, _ = curve_fit_wrapper(f, time_s, m, (1, 1, 1), bounds=bounds)
        else:
            raise e
    return popt


def timetable(timetable_path: Path):
    timetable = pd.read_csv(timetable_path)
    t = np.array(
        [
            timetable.loc[
                (timetable["sequence_label"] == "mixed")
                & (timetable["subject"] == "sub-01")
                & (timetable["session"] == session)
            ]["acquisition_relative_injection"].item()
            for session in [f"ses-{i:02d}" for i in range(1, 6)]
        ]
    )
    return np.maximum(0, t) / (3600 * 24)


def fit_voxel_wrapper(t: np.ndarray, pbar, f):
    def call(y):
        tj = 1.2 * t[np.argmax(y)]
        b = 1
        c = tj / b
        a = 1.2 * y.max() / (f(b * c, 1, b, c))
        p0 = (a, b, c)
        return fit_voxel(t, pbar, y, f, p0, bounds=([0, 0.1, 0], [100, 20, 20]))

    return call


def curve_fit_csf_concentrations(
    concentrations: list[Path],
    csf_mask_path: Path,
    timetable_path: Path,
) -> list[nibabel.nifti1.Nifti1Image]:
    D = np.array(
        [nibabel.nifti1.load(ci).get_fdata(dtype=np.single) for ci in concentrations]
    )
    csf_mask_nii = nibabel.nifti1.load(csf_mask_path)
    csf_mask = csf_mask_nii.get_fdata(dtype=np.single).astype(bool)
    valid_voxels = np.array(np.where(csf_mask)).T
    Dmasked = np.array([D[:, i, j, k] for (i, j, k) in valid_voxels])

    t = timetable(timetable_path)
    with tqdm.tqdm(total=len(Dmasked)) as pbar:
        voxel_fitter = fit_voxel_wrapper(t, pbar, f)
        vfunc = np.vectorize(voxel_fitter, signature="(n) -> (3)")
        fitted_coefficients = vfunc(Dmasked)

    a, b, c = (
        fitted_coefficients[:, 0],
        fitted_coefficients[:, 1],
        fitted_coefficients[:, 2],
    )

    I, J, K = valid_voxels.T  # noqa: E741
    coefficient_map = np.nan * np.zeros((3, *csf_mask.shape))
    coefficient_map[0, I, J, K] = a
    coefficient_map[1, I, J, K] = b
    coefficient_map[2, I, J, K] = c

    affine = nibabel.nifti1.load(csf_mask_path).affine
    return [nibabel.nifti1.Nifti1Image(im, affine) for im in coefficient_map]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--inputfiles", nargs="+", type=Path, required=True)
    parser.add_argument("--csfmask", type=Path, required=True)
    parser.add_argument("--timetable", type=Path, required=True)
    parser.add_argument("--a_out", type=Path, required=True)
    parser.add_argument("--b_out", type=Path, required=True)
    parser.add_argument("--c_out", type=Path, required=True)
    args = parser.parse_args()

    A, B, C = curve_fit_csf_concentrations(
        args.inputfiles, args.csfmask, args.timetable
    )
    for im, path in zip((A, B, C), (args.a_out, args.b_out, args.c_out)):
        nibabel.nifti1.save(im, path)
