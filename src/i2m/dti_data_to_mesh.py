"""
Adjusted version of script from original MRI2FEM-book according to below license.
'''
Copyright (c) 2020 Kent-Andre Mardal, Marie E. Rognes, Travis B. Thompson, Lars Magnus Valnes

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''
"""

from pathlib import Path

import click
import dolfin as df
import nibabel
import numpy as np
import panta_rhei as pr
import simple_mri as sm
from nibabel.affines import apply_affine
from ufl import dev, inner, sqrt, tr

from dti.clean_dti_data import extend_to_9_component_array
from i2m.mri2fenics import find_dof_nearest_neighbours, locate_dof_voxels


def mean_diffusivity(Dvector: np.ndarray) -> np.ndarray:
    return Dvector[..., ::4].sum() / 3.0


def adjusting_mean_diffusivity(Dvector: np.ndarray, subdomains, tags_with_limits):
    MD = mean_diffusivity(Dvector)
    # Reads the tag and the minimum and maximum mean diffusivity limit
    # for that tag subdomain.
    for tag, mn, mx in tags_with_limits:
        # If the minimum or maximum mean diffusivity limit is set to zero,
        # then the limit is considered void.
        usr_max = float(mx) if mx != 0 else np.inf
        usr_min = float(mn) if mn != 0 else -np.inf

        # creates a mask for all degrees of freesom that are within the
        # subdomain with the tag and is above the maximum limit or
        # below the minimum limit.
        max_mask = (subdomains.array() == tag) * (MD > usr_max)
        min_mask = (subdomains.array() == tag) * (MD < usr_min)

        # Sets values that are either above or below limits to the closest limit.
        Dvector[max_mask] = usr_max * np.divide(
            Dvector[max_mask], MD[max_mask, np.newaxis]
        )
        Dvector[min_mask] = usr_min * np.divide(
            Dvector[min_mask], MD[min_mask, np.newaxis]
        )


def dti_data_to_mesh(
    domain: pr.Domain,
    dti_mri: sm.SimpleMRI,
    md_mri: sm.SimpleMRI,
    fa_mri: sm.SimpleMRI,
    brain_mask: np.ndarray,
    output: Path,
):
    # Structure tensors as 9-component row-major vectors
    dti = dti_mri.data.squeeze()
    dti = extend_to_9_component_array(dti)
    valid_dti = np.isfinite(dti).all(axis=-1)
    dti_mask = valid_dti * brain_mask

    DG0 = df.FunctionSpace(domain, "DG", 0)
    DG09 = df.TensorFunctionSpace(domain, "DG", 0)

    dof_voxels = locate_dof_voxels(DG0, md_mri, rint=False)
    dof_neighbours = find_dof_nearest_neighbours(dof_voxels, dti_mask, N=1)

    D = df.Function(DG09)
    D.vector()[:] = dti[*dof_neighbours].ravel()

    md_mask = (md_mri.data > 1e-5) * brain_mask
    dof_N_neighbours = find_dof_nearest_neighbours(dof_voxels, md_mask, N=10)
    MD = df.Function(DG0)
    MD.vector()[:] = np.median(md_mri.data[*dof_N_neighbours], axis=0)

    FA = df.Function(DG0)
    FA.vector()[:] = np.median(fa_mri.data[*dof_N_neighbours], axis=0)

    hdf = df.HDF5File(domain.mpi_comm(), str(output), "w")
    pr.write_domain(hdf, domain)
    pr.write_function(hdf, D, "DTI")
    pr.write_function(hdf, MD, "MD")
    pr.write_function(hdf, FA, "FA")
    hdf.close()


@click.command()
@click.option("--mesh", type=Path, required=True)
@click.option("--dti", type=Path, required=True)
@click.option("--md", type=Path, required=True)
@click.option("--fa", type=Path, required=True)
@click.option("--mask", type=Path, required=True)
@click.option("--output", type=Path, required=True)
def dti2mesh(mesh: Path, dti: Path, md: Path, fa: Path, mask: Path, output: Path):
    dti_mri = sm.load_mri(dti, dtype=np.double)
    md_mri = sm.load_mri(md, dtype=np.double)
    fa_mri = sm.load_mri(fa, dtype=np.double)
    mask_mri = sm.load_mri(mask, dtype=bool)

    hdf = df.HDF5File(df.MPI.comm_world, str(mesh), "r")
    domain = pr.read_domain(hdf)
    hdf.close()
    dti_data_to_mesh(domain, dti_mri, md_mri, fa_mri, mask_mri.data, output)


if __name__ == "__main__":
    dti2mesh()
