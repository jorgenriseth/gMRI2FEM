import re
from pathlib import Path
from typing import Literal, Optional

import numpy as np
import pandas as pd
import scipy as sp
import skimage
import tqdm


def apply_affine(T: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Apply homogeneous-coordinate affine matrix T to each row of of matrix
    X of shape (N, 3)"""
    A = T[:-1, :-1]
    b = T[:-1, -1]
    return A.dot(X.T).T + b


def threshold_between(
    x: float | np.ndarray, lo: float, hi: float
) -> float | np.ndarray:
    return np.maximum(lo, np.minimum(x, hi))


def nan_filter_gaussian(
    U: np.ndarray, sigma: float, truncate: float = 4.0
) -> np.ndarray:
    V = U.copy()
    V[np.isnan(U)] = 0
    VV = sp.ndimage.gaussian_filter(V, sigma=sigma, truncate=truncate)

    W = np.ones_like(U)
    W[np.isnan(U)] = 0
    WW = sp.ndimage.gaussian_filter(W, sigma=sigma, truncate=truncate)
    mask = ~((WW == 0) * (VV == 0))
    out = np.nan * np.zeros_like(U)
    out[mask] = VV[mask] / WW[mask]
    return out


def smooth_extension(D: np.ndarray, sigma: float, truncate: float = 4) -> np.ndarray:
    return np.where(np.isnan(D), nan_filter_gaussian(D, sigma, truncate), D)


def mri_facemask(vol: np.ndarray, smoothing_level=5):
    thresh = skimage.filters.threshold_triangle(vol)
    binary = vol > thresh
    binary = sp.ndimage.binary_fill_holes(binary)
    binary = skimage.filters.gaussian(binary, sigma=smoothing_level)
    binary = binary > skimage.filters.threshold_isodata(binary)
    return binary


def largest_island(mask: np.ndarray, connectivity: int = 1) -> np.ndarray:
    newmask = skimage.measure.label(mask, connectivity=connectivity)
    regions = skimage.measure.regionprops(newmask)
    regions.sort(key=lambda x: x.num_pixels, reverse=True)
    return newmask == regions[0].label


def grow_restricted(grow, restriction, growth_radius):
    return (
        grow
        + skimage.morphology.binary_dilation(
            grow, skimage.morphology.cube(2 * growth_radius + 1)
        )
        * restriction
    )


def segmentation_smoothing(
    segmentation, sigma, cutoff_score=0.5, **kwargs
) -> dict[str, np.ndarray]:
    labels = np.unique(segmentation)
    labels = labels[labels != 0]
    new_labels = np.zeros_like(segmentation)
    high_scores = np.zeros(segmentation.shape)
    for label in tqdm.tqdm(labels):
        label_scores = sp.ndimage.gaussian_filter(
            (segmentation == label).astype(float), sigma=sigma, **kwargs
        )
        is_new_high_score = label_scores > high_scores
        new_labels[is_new_high_score] = label
        high_scores[is_new_high_score] = label_scores[is_new_high_score]

    delete_scores = (high_scores < cutoff_score) * (segmentation == 0)
    new_labels[delete_scores] = 0
    return {"labels": new_labels, "scores": high_scores}


def nearest_neighbour(
    D: np.ndarray, inds: np.ndarray, valid_indices: Optional[np.ndarray] = None
) -> np.ndarray:
    i, j, k = inds.T
    if valid_indices is None:
        I, J, K = np.array(np.where(np.isfinite(D)))  # noqa: E741
    else:
        I, J, K = valid_indices  # noqa: E741
    interp = sp.interpolate.NearestNDInterpolator(np.array((I, J, K)).T, D[I, J, K])
    D_out = D.copy()
    D_out[i, j, k] = interp(i, j, k)
    num_nan_values = (~np.isfinite(D_out[i, j, k])).sum()
    assert num_nan_values == 0
    return D_out


def plot_orient(
    volume: np.ndarray, plane: Literal["sagittal", "coronal", "axial"], idx: int
) -> np.ndarray:
    if plane == "sagittal":
        return volume[idx, :, ::-1].T
    elif plane == "coronal":
        return volume[:, idx, ::-1].T
    elif plane == "axial":
        return volume[:, ::-1, idx].T
    raise ValueError(
        f"Invalid plane '{plane}', should be one of 'sagittal', 'coronal', 'axial'"
    )


def session_range(*args):
    return [f"ses-{idx + 1:02d}" for idx in range(*args)]


def axial_dilation_footprint(ndim, axis):
    footprint = [1] * ndim
    footprint[axis] = 3
    return np.ones(footprint)


def connectivity_matrix(arr):
    labels = np.unique(arr)
    n_labels = len(labels)
    K = np.zeros((n_labels, n_labels))
    for i, label_i in enumerate(labels):
        mask_i = arr == label_i
        axially_dilated_masks = [
            skimage.morphology.dilation(
                mask_i, footprint=np.array(axial_dilation_footprint(arr.ndim, axis))
            )
            for axis in range(arr.ndim)
        ]
        K[i, i] = mask_i.sum()
        for j, label_j in enumerate(labels[i + 1 :], start=i + 1):
            mask_j = arr == label_j
            interfaces = [(dilated * mask_j).sum() for dilated in axially_dilated_masks]
            print(interfaces)
            K[i, j] = sum(interfaces)

    K += np.triu(K, k=1).T
    return {"labels": labels, "connectivity": K}


def float_string_formatter(x: float, digits):
    if float(x) == float("inf"):
        return "inf"
    return f"{x * 10 ** (-digits):{f'.{digits}e'}}".replace(".", "")


def to_scientific(num, decimals):
    if float(num) == float("inf"):
        return r"\infty"
    x = f"{float(num):{f'.{decimals}e'}}"
    m = re.search(r"(\d\.{0,1}\d*)e([\+|\-]\d{2})", x)

    return f"{m.group(1)}\\times10^{{{int(m.group(2))}}}"


def find_timestamps(
    timetable_path: Path,
    sequence_name: str,
    subject: str,
):
    try:
        timetable = pd.read_csv(timetable_path, sep="\t")
    except pd.errors.EmptyDataError:
        raise RuntimeError(f"Timetable-file {timetable_path} is empty.")
    if "sequence_name" in timetable.columns:
        seqlabel = "sequence_name"
    elif "sequence_label" in timetable.columns:
        seqlabel = "sequence_label"
    else:
        raise RuntimeError("Cant find column 'sequence_name' or 'sequence_label'")
    subject_sequence_entries = (
        (timetable.subject == subject) 
        &(timetable[seqlabel].str.lower() == sequence_name)
    )  # fmt: skip
    try:
        acq_times = timetable.loc[subject_sequence_entries][
            "acquisition_relative_injection"
        ]
        times = np.array(acq_times)
        assert len(times) > 0, f"Couldn't find time for {subject}: {sequence_name}"
    except (AssertionError, ValueError) as e:
        print(timetable)
        print(subject, sequence_name)
        print(subject_sequence_entries)
        raise e
    return times


def find_timestamp(
    timetable_path: Path,
    timestamp_sequence: str,
    subject: str,
    session: str,
) -> float:
    """Find single session timestamp"""
    try:
        timetable = pd.read_csv(timetable_path, sep="\t")
    except pd.errors.EmptyDataError:
        raise RuntimeError(f"Timetable-file {timetable_path} is empty.")
    try:
        timestamp = timetable.loc[
            (timetable["sequence_name"].str.lower() == timestamp_sequence)
            & (timetable["subject"] == subject)
            & (timetable["session"] == session)
        ]["acquisition_relative_injection"]
    except ValueError as e:
        print(timetable)
        print(timestamp_sequence, subject)
        raise e
    return timestamp.item()


def prepend_info(df, **kwargs):
    nargs = len(kwargs)
    for key, val in kwargs.items():
        assert key not in df.columns, f"Column {key} already exist in df."
        df[key] = val
    return df[[*df.columns[-nargs:], *df.columns[:-nargs]]]


def with_suffix(p: Path, newsuffix: str) -> Path:
    return p.parent / f"{p.name.split('.')[0]}{newsuffix}"
