from pathlib import Path

import click
import dolfin as df
import numpy as np
import simple_mri as sm
import tqdm
from pantarei import read_domain
from simple_mri import load_mri

from gmri2fem.utils import apply_affine
from i2m.concentrations_to_mesh import nearest_neighbour
from i2m.mri2fenics import find_dof_nearest_neighbours, locate_dof_voxels


@click.command("subdomains")
@click.argument("segmentation", type=Path, required=True)
@click.argument("meshpath", type=Path)
@click.argument("output", type=Path)
def subdomains(*args, **kwargs):
    segments_to_mesh(*args, **kwargs)


CEREBRAL_WM_RANGES = [
    *[2, 41],  # aseg left/right cerebral white labels
    *list(range(3000, 3036)),  # wmparc-left-labels
    *list(range(4000, 4036)),  # wmparc-right-labels
    *[5000, 5001],
    *[28, 60],  # VentralDC included in white matter sudbomain
    *list(range(251, 256)),  #  Corpus callosum
    *[31, 63],  # Choroid plexus.
]
CEREBRAL_GM_RANGES = [
    *[3, 42],  # aseg left/right cortcal gm
    *list(range(1000, 1036)),  # aparc left labels
    *list(range(2000, 2036)),  # aparc right labels
]
SUBCORTICAL_GM_RANGES = [
    *(10, 49),  # Thalamus,
    *(11, 50),  # Caudate,
    *(12, 51),  # Putamen,
    *(13, 52),  # pallidum
    *(17, 53),  # hippocampus
    *(18, 54),  # amygdala
    *(26, 58),  # accumbens
]
subdomain_allowed = {
    1: CEREBRAL_GM_RANGES,
    2: CEREBRAL_WM_RANGES,
    3: SUBCORTICAL_GM_RANGES,
}


def extract_label_modes(seg_NN):
    return np.apply_along_axis(lambda x: np.bincount(x).argmax(), axis=0, arr=seg_NN)


def segments_to_mesh(
    segmentation: Path,
    meshpath: Path,
    output: Path,
):
    hdf = df.HDF5File(df.MPI.comm_world, str(meshpath), "r")
    domain = read_domain(hdf)
    mesh, subdomains, _ = domain, domain.subdomains, domain.boundaries
    d = mesh.topology().dim()

    segmentation_mri = load_mri(segmentation, dtype=int)
    seg = segmentation_mri.data

    n = 3**d
    N = mesh.num_cells()

    cell_midpoints = np.array([cell.midpoint()[:] for cell in df.cells(mesh)])
    cell_indices = np.array([cell.index() for cell in df.cells(mesh)])
    regions = np.zeros(N)
    for subdomain_label in np.unique(subdomains.array()):
        subdomain_mask = np.isin(seg, subdomain_allowed[subdomain_label])
        subdomain_cells = subdomains.array() == subdomain_label
        subdomain_cell_midpoints = cell_midpoints[subdomain_cells]
        subdomain_cell_NN = find_dof_nearest_neighbours(
            sm.apply_affine(
                np.linalg.inv(segmentation_mri.affine), subdomain_cell_midpoints
            ),
            subdomain_mask,
            n,
        )
        seg_NN = seg[*subdomain_cell_NN]
        regions[subdomain_cells] = extract_label_modes(seg_NN)

    parcellations = df.MeshFunction("size_t", mesh, d, 0)
    parcellations.array()[cell_indices] = regions[cell_indices]

    hdf = df.HDF5File(mesh.mpi_comm(), str(output.with_suffix(".hdf")), "w")
    hdf.write(parcellations, "/parcellations")
    hdf.close()


if __name__ == "__main__":
    subdomains()
