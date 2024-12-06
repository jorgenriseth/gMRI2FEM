import re
from pathlib import Path
from typing import Optional


import click
import dolfin as df
import numpy as np
import scipy
import pandas as pd
import tqdm
import pantarei as pr
import skimage
from pantarei import FenicsStorage, fenicsstorage2xdmf

from gmri2fem.utils import apply_affine, nan_filter_gaussian
from simple_mri import load_mri, SimpleMRI
from i2m.mri2fenics import (
    locate_dof_voxels,
    find_boundary_dofs,
    find_dof_nearest_neighbours,
)


def extract_sequence_timestamps(
    timetable_path: Path, subject: str, sequence_label: str
):
    try:
        timetable = pd.read_csv(timetable_path, sep="\t")
    except pd.errors.EmptyDataError:
        raise RuntimeError(f"Timetable-file {timetable_path} is empty.")
    subject_sequence_entries = (timetable["subject"] == subject) & (
        timetable["sequence_label"] == sequence_label
    )
    try:
        acq_times = timetable[subject_sequence_entries][
            "acquisition_relative_injection"
        ]
    except ValueError as e:
        print(timetable)
        print(subject, sequence_label)
        print(subject_sequence_entries)
        raise e
    times = np.array(acq_times)
    assert len(times) > 0, f"Couldn't find time for {subject}: {sequence_label}"
    return times


def smooth_dilation(D: np.ndarray, sigma: float, truncate: float = 4) -> np.ndarray:
    return np.where(np.isfinite(D), D, nan_filter_gaussian(D, sigma, truncate))


def nearest_neighbour(
    D: np.ndarray, inds: np.ndarray, valid_indices: Optional[np.ndarray] = None
) -> np.ndarray:
    i, j, k = inds.T
    if valid_indices is None:
        I, J, K = np.array(np.where(np.isfinite(D)))
    else:
        I, J, K = valid_indices
    interp = scipy.interpolate.NearestNDInterpolator(np.array((I, J, K)).T, D[I, J, K])
    D_out = D.copy()
    D_out[i, j, k] = interp(i, j, k)
    num_nan_values = (~np.isfinite(D_out[i, j, k])).sum()
    assert num_nan_values == 0
    return D_out


def smooth_extension(D, inds, sigma, truncate, maxiter=10, fallback=True):
    i, j, k = inds.T
    num_nan_values = (~np.isfinite(D[i, j, k])).sum()
    for _ in range(maxiter):
        D = smooth_dilation(D, sigma, truncate)
        num_nan_values = (~np.isfinite(D[i, j, k])).sum()
        if num_nan_values == 0:
            return D

    if not fallback:
        raise RuntimeError(
            f"Couldn't extend data within {maxiter} iterations, missing {num_nan_values}"
        )
    return nearest_neighbour(D, inds)


def map_concentration(
    subject: str,
    concentration_paths: list[Path],
    meshpath: Path,
    csfmask_path: Path,
    timetable: Path,
    output: Path,
    femfamily: str,
    femdegree: int,
    visualdir: Optional[Path] = None,
):
    csf_mask_mri = load_mri(csfmask_path, dtype=bool)
    csf_mask = skimage.morphology.binary_erosion(
        csf_mask_mri.data, skimage.morphology.ball(1)
    )

    timestamps = np.maximum(
        0, extract_sequence_timestamps(timetable, subject, "looklocker")
    )
    domain = pr.hdf2fenics(meshpath, pack=True)
    V = df.FunctionSpace(domain, femfamily, femdegree)

    concentration_mri = load_mri(concentration_paths[0], dtype=np.single)
    dof_voxels = locate_dof_voxels(V, concentration_mri, rint=False)

    boundary_dofs = find_boundary_dofs(V)
    boundary_dof_neighbours = find_dof_nearest_neighbours(
        dof_voxels[boundary_dofs], csf_mask, N=10
    )
    assert len(concentration_paths) > 0
    assert len(timestamps) == len(concentration_paths)

    outfile = FenicsStorage(str(output), "w")
    outfile.write_domain(domain)
    for ti, ci in zip(tqdm.tqdm(timestamps), concentration_paths):
        concentration_mri = load_mri(ci, dtype=np.single)
        valid_concentrations = np.isfinite(concentration_mri.data)
        boundary_dof_neighbours = find_dof_nearest_neighbours(
            dof_voxels[boundary_dofs],
            csf_mask * valid_concentrations,
            N=10,
        )
        u_boundary = df.Function(V)
        u_boundary.vector()[boundary_dofs] = np.nanmedian(
            concentration_mri.data[*boundary_dof_neighbours], axis=0
        )
        outfile.write_checkpoint(u_boundary, name="boundary_concentration", t=ti)

        dof_neighbours = find_dof_nearest_neighbours(
            dof_voxels, ~csf_mask * valid_concentrations, N=10
        )
        u_internal = df.Function(V)
        u_internal.vector()[:] = np.nanmedian(
            concentration_mri.data[*dof_neighbours], axis=0
        )
        outfile.write_checkpoint(u_internal, name="concentration", t=ti)
    outfile.close()

    if visualdir is not None:
        fenicsstorage2xdmf(
            FenicsStorage(outfile.filepath, "r"),
            "concentration",
            "internal",
            lambda _: visualdir / "concentrations_internal.xdmf",
        )
        fenicsstorage2xdmf(
            FenicsStorage(outfile.filepath, "r"),
            "boundary_concentration",
            "boundary",
            lambda _: visualdir / "concentrations_boundary.xdmf",
        )


@click.command()
@click.argument("concentration_paths", type=Path, nargs=-1, required=True)
@click.option("--meshpath", type=Path, required=True)
@click.option("--csfmask_path", type=Path, required=True)
@click.option("--timetable", type=Path, required=True)
@click.option("--output", type=Path, required=True)
@click.option("--femfamily", type=str, default="CG")
@click.option("--femdegree", type=int, default=1)
@click.option("--visualdir", type=Path)
def concentrations2mesh(concentration_paths, **kwargs):
    subject_re = re.compile(r"(?P<subject>sub-(control|patient)*\d{2})")
    m = subject_re.search(str(concentration_paths[0]))
    if m is None:
        raise ValueError(f"Couldn't find subject in path {concentration_paths[0]}")
    map_concentration(m.groupdict()["subject"], concentration_paths, **kwargs)


if __name__ == "__main__":
    concentrations2mesh()
