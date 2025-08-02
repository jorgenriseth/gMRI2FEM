import re
from pathlib import Path
from typing import Optional

import click
import dolfin as df
import numpy as np
import panta_rhei as pr
import skimage
import tqdm
import simple_mri as sm
from panta_rhei import FenicsStorage, fenicsstorage2xdmf
from simple_mri import load_mri
from loguru import logger

from gmri2fem.utils import find_timestamps, nan_filter_gaussian, nearest_neighbour
from i2m.mri2fenics import (
    find_boundary_dofs,
    find_dof_nearest_neighbours,
    locate_dof_voxels,
    mri2fem_interpolate_quadrature,
)


def smooth_dilation(D: np.ndarray, sigma: float, truncate: float = 4) -> np.ndarray:
    return np.where(np.isfinite(D), D, nan_filter_gaussian(D, sigma, truncate))


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
    collocation: bool = False,
    collocation_npoints: int = 10,
    quad_degree: int = 6,
):
    csf_mask_mri = load_mri(csfmask_path, dtype=bool)
    csf_mask = skimage.morphology.binary_erosion(
        csf_mask_mri.data, skimage.morphology.ball(1)
    )

    timestamps = np.maximum(0, find_timestamps(timetable, "looklocker", subject))
    domain = pr.hdf2fenics(meshpath, pack=True)
    V = df.FunctionSpace(domain, femfamily, femdegree)

    concentration_mri = load_mri(concentration_paths[0], dtype=np.single)
    dof_voxels = locate_dof_voxels(V, concentration_mri, rint=False)

    boundary_dofs = find_boundary_dofs(V)
    assert len(concentration_paths) > 0
    assert len(timestamps) == len(concentration_paths)

    outfile = FenicsStorage(str(output), "w")
    outfile.write_domain(domain)
    for ti, ci in zip(tqdm.tqdm(timestamps), concentration_paths):
        concentration_mri = load_mri(ci, dtype=np.single)
        valid_concentrations = np.isfinite(concentration_mri.data)
        u_boundary = map_boundary_concentration(
            concentration_mri,
            csf_mask * valid_concentrations,
            V,
            boundary_dofs,
            dof_voxels,
            collocation_npoints,
        )
        outfile.write_checkpoint(u_boundary, name="boundary_concentration", t=ti)

        if collocation:
            dof_neighbours = find_dof_nearest_neighbours(
                dof_voxels, ~csf_mask * valid_concentrations, N=collocation_npoints
            )
            u_internal = df.Function(V)
            u_internal.vector()[:] = np.nanmedian(
                concentration_mri.data[*dof_neighbours], axis=0
            )
        else:
            try:
                u_internal = mri2fem_interpolate_quadrature(
                    concentration_mri, V, quad_degree, mask=(~csf_mask)
                )
            except RuntimeError:
                logger.warning(
                    "Galerkin projection failed, probably due to memory. Falling back to CG-solver"
                )
                u_internal = mri2fem_interpolate_quadrature(
                    concentration_mri,
                    V,
                    quad_degree,
                    mask=(~csf_mask),
                    solver_type="iterative",
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


@click.command("boundary-concentrations")
@click.argument("mris", type=str, nargs=-1)
@click.option("--output", "-o", type=str, required=True)
@click.option("--mesh", "-m", type=str, required=True)
@click.option("--csfmask", "-c", type=str, required=True)
@click.option("--timestamps", "-t", type=str, required=True)
@click.option("--npoints", type=int, default=27)
def map_boundary_concentrations_cli(
    mris, mesh: str, csfmask: str, timestamps: str, output: str, npoints: int = 27
):
    time = np.loadtxt(timestamps)
    assert len(mris) > 0, "No MRI's provided"
    assert len(time) == len(mris), "Timestamps and MRI's have different length"

    csf_mask_mri = load_mri(csfmask, dtype=bool)
    csf_mask = skimage.morphology.binary_erosion(
        csf_mask_mri.data, skimage.morphology.ball(1)
    )
    domain = pr.hdf2fenics(mesh, pack=True)
    V = df.FunctionSpace(domain, "CG", 1)

    reference_mri = load_mri(mris[0], dtype=np.single)
    dof_voxels = locate_dof_voxels(V, reference_mri, rint=False)

    boundary_dofs = find_boundary_dofs(V)
    boundary_functions = [
        df.Function(V, name="concentration") for _ in range(len(mris))
    ]
    for idx, ci in enumerate(tqdm.tqdm(mris)):
        concentration_mri = load_mri(ci, dtype=np.single)
        valid_concentrations = np.isfinite(concentration_mri.data)
        boundary_functions[idx].assign(
            map_boundary_concentration(
                concentration_mri,
                csf_mask * valid_concentrations,
                V,
                boundary_dofs,
                dof_voxels,
                npoints,
            )
        )
    with df.HDF5File(domain.mpi_comm(), str(output), "w") as hdf:
        pr.write_domain(hdf, domain)
        pr.write_function(hdf, boundary_functions[0], "boundary_concentration")
        for ti, ui in zip(time[1:], boundary_functions[1:]):
            pr.write_checkpoint(hdf, ui, name="boundary_concentration", t=float(ti))


def map_boundary_concentration(
    concentration_mri: sm.SimpleMRI,
    csf_mask: np.ndarray,
    V: df.FunctionSpace,
    boundary_dofs: np.ndarray,
    dof_voxels: np.ndarray,
    npoints: int,
):
    boundary_dof_neighbours = find_dof_nearest_neighbours(
        dof_voxels[boundary_dofs],
        csf_mask,
        N=npoints,
    )
    u_boundary = df.Function(V)
    u_boundary.vector()[boundary_dofs] = np.nanmedian(
        concentration_mri.data[*boundary_dof_neighbours], axis=0
    )
    return u_boundary


@click.command()
@click.argument("concentration_paths", type=Path, nargs=-1, required=True)
@click.option("--meshpath", type=Path, required=True)
@click.option("--csfmask_path", type=Path, required=True)
@click.option("--timetable", type=Path, required=True)
@click.option("--output", type=Path, required=True)
@click.option("--femfamily", type=str, default="CG")
@click.option("--femdegree", type=int, default=1)
@click.option("--subject_regex", type=str)
@click.option("--visualdir", type=Path)
def concentrations2mesh(concentration_paths, subject_regex, **kwargs):
    subject_regex = subject_regex or r"sub-(control|patient)*\d{2}"
    subject_re = re.compile(rf"(?P<subject>{subject_regex})")
    m = subject_re.search(str(concentration_paths[0]))
    if m is None:
        raise ValueError(f"Couldn't find subject in path {concentration_paths[0]}")
    map_concentration(
        m.groupdict()["subject"], concentration_paths, **kwargs, collocation=False
    )


if __name__ == "__main__":
    concentrations2mesh()
