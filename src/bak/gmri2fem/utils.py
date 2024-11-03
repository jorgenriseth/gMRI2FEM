import json
import re
import shutil
from datetime import datetime
from functools import partial
from itertools import chain, repeat
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Iterator
import skimage

import numpy as np
import pydicom
import scipy as sp


def create_protocol_filemap(
    input_dir: Path, output_dir: Path, sequences: dict[str, str], n_jobs=None
) -> dict[Path, Path]:
    """Runs through the folder structure with the MRI data we receive from the hospital, and creates a dictionary mapping a src-path sorts them according to
    'patient - study_datetime - protocol'. Runs in parallel for each of the studies."""
    pool = Pool(n_jobs)
    task = partial(
        create_study_filemap,
        output_dir=output_dir,
        sequences=sequences,
    )
    results = pool.map(task, study_iterator(input_dir))
    return dict(chain(*(x.items() for x in results)))


def create_study_filemap(
    study_dir: Path, output_dir: Path, sequences: dict[str, str]
) -> dict[Path, Path]:
    date_dir = study_dir.parent
    patient = date_dir.parent.stem
    study_data = study_dir / "DICOM" / "DICOM"
    date = datetime.strptime(date_dir.stem, "%Y_%m_%d").strftime("%Y%m%d")
    timestamp = find_timestamp(study_dir)
    study_target = output_dir / patient / f"{date}_{timestamp}"

    filemap = {}
    for imfile, offset in study_imfiles(study_data):
        with pydicom.dcmread(imfile) as f:
            file_protocol = f.ProtocolName
            if file_protocol in sequences:
                filemap[imfile] = (
                    study_target
                    / sequences[file_protocol]
                    / renumber_imfile(imfile, offset)
                )
    return filemap


def find_timestamp(study_dir: Path) -> str:
    study_data_path = study_dir / "DICOM" / "DICOM"
    first_imfile, _ = next(study_imfiles(study_data_path))
    with pydicom.dcmread(first_imfile) as f:
        timestamp = f.StudyTime
    return timestamp


def study_imfiles(study_data_path: Path) -> Iterator[Path]:
    return chain(
        study_imfiles_in_dicomdir(study_data_path),
        *(
            study_imfiles_in_dicom_subdir(subdir)
            for subdir in filter(is_image_subdirectory, study_data_path.iterdir())
        ),
    )


def study_imfiles_in_dicomdir(study_data_path: Path) -> Iterator[tuple[Path, int]]:
    return zip(filter(is_imfile, study_data_path.iterdir()), repeat(0))


def study_imfiles_in_dicom_subdir(
    study_subdir_path: Path,
) -> Iterator[tuple[Path, int]]:
    return zip(
        filter(is_imfile, study_subdir_path.iterdir()),
        repeat(int(study_subdir_path.stem)),
    )


def study_iterator(input_dir: Path) -> Iterator[Path]:
    return (
        x
        for patient in input_dir.iterdir()
        for date in patient.iterdir()
        for x in date.iterdir()
    )


def renumber_imfile(imfile: Path, offset: int) -> str:
    return f"IM_{int(imfile.stem.split('_')[1]) + offset * 2048:04d}"


def is_image_subdirectory(path: Path) -> bool:
    return path.is_dir() and path.stem.isdigit()


def is_imfile(path: Path) -> bool:
    """Checks if path is a dicom IM-file."""
    return path.is_file() and path.stem[:3] == "IM_"


def create_study_metadata(study_path: Path):
    study_metadata = {}
    for path in filter(lambda x: x.is_dir(), study_path.iterdir()):
        protocol_filelist = [p for p in path.iterdir() if is_imfile(p)]
        label_list = [int(x.stem.split("_")[1]) for x in protocol_filelist]
        with pydicom.dcmread(protocol_filelist[0]) as f:
            timestamp = f.SeriesTime

        study_metadata[path.stem] = {
            "min": min(label_list),
            "max": max(label_list),
            "num_files": len(protocol_filelist),
            "series_time": timestamp,
        }
    return study_metadata


def store_study_metadata(study_path: Path) -> None:
    filepath = study_path / "info.json"
    with open(filepath, "w") as f:
        json.dump(create_study_metadata(study_path), f, indent=4)


def add_sorted_metadata(output_dir: Path) -> None:
    reorganized_study_iterator = (
        study for patient in output_dir.iterdir() for study in patient.iterdir()
    )
    for study in reorganized_study_iterator:
        store_study_metadata(study)


def float_string_formatter(x: float, digits):
    if float(x) == float("inf"):
        return "inf"
    return f"{x*10**(-digits):{f'.{digits}e'}}".replace(".", "")


def to_scientific(num, decimals):
    if float(num) == float("inf"):
        return "\infty"
    x = f"{float(num):{f'.{decimals}e'}}"
    m = re.search("(\d\.{0,1}\d*)e([\+|\-]\d{2})", x)

    return f"{m.group(1)}\\times10^{{{int(m.group(2))}}}"


def nested_dict_set(d: dict[str, dict | float], keys: tuple[str], value: float):
    if isinstance(keys, str):
        d[keys] = value
        return d

    depth = len(keys)
    d_ = d
    for i, key in enumerate(keys):
        if i == depth - 1:
            d_[key] = value
        else:
            d_ = d[key]
    return d


def apply_affine(T: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Apply homogeneous-coordinate affine matrix T to each row of of matrix
    X of shape (N, 3)"""
    A = T[:-1, :-1]
    b = T[:-1, -1]
    return A.dot(X.T).T + b


def threshold_between(x: float | np.ndarray, lo: float, hi: float):
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
