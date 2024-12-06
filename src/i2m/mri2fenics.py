from pathlib import Path
from typing import Callable, Optional

import dolfin as df
import nibabel
import numpy as np
import scipy
import simple_mri as sm
from pantarei import FenicsStorage

from gmri2fem.utils import apply_affine, nan_filter_gaussian


def find_dof_nearest_neighbours(
    dof_inds: np.ndarray, mask: np.ndarray, N: int
) -> np.ndarray:
    dof_neighbours = -np.ones((2, N, dof_inds.shape[0]), int)
    valid_inds = np.argwhere(mask)
    tree = scipy.spatial.KDTree(valid_inds)
    distances, indices = tree.query(dof_inds, k=N)
    dof_neighbours = valid_inds[indices].T
    return dof_neighbours


def find_boundary_dofs(V: df.FunctionSpace) -> np.ndarray:
    return np.array(
        [
            dof
            for dof in df.DirichletBC(V, df.Constant(0), "on_boundary")
            .get_boundary_values()
            .keys()
        ]
    )


def locate_dof_voxels(V: df.FunctionSpace, mri: sm.SimpleMRI, rint: bool = True):
    """Create a list of indices of voxels of an mri for which the dof coordinates
    of a fenics function space are located within."""
    dof_coordinates = V.tabulate_dof_coordinates()
    img_space_coords = sm.apply_affine(np.linalg.inv(mri.affine), dof_coordinates)
    if rint:
        return np.rint(img_space_coords).astype(int)
    return img_space_coords


def mri2fem_interpolate(
    D: np.ndarray,
    affine: np.ndarray,
    V: df.FunctionSpace,
    datafilter: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> df.Function:
    if datafilter is not None:
        D = datafilter(D)
    u = df.Function(V)
    z = V.tabulate_dof_coordinates()
    ind = np.rint(apply_affine(np.linalg.inv(affine), z)).astype(int)
    i, j, k = ind.T
    u.vector()[:] = D[i, j, k]
    return u


def smooth_extension(D: np.ndarray, sigma: float, truncate: float = 4) -> np.ndarray:
    return np.where(np.isnan(D), nan_filter_gaussian(D, sigma, truncate), D)


def read_image(
    filename: Path,
    functionspace: df.FunctionSpace,
    datafilter: Optional[Callable[[np.ndarray], np.ndarray]] = None,
) -> df.Function:
    mri_volume = nibabel.nifti1.load(filename)
    voxeldata = mri_volume.get_fdata(dtype=np.single)
    return mri2fem_interpolate(voxeldata, mri_volume.affine, functionspace, datafilter)


def fenicsstorage2xdmf(
    filepath, funcname: str, subnames: str | list[str], outputdir: Path
) -> None:
    file = FenicsStorage(filepath, "r")
    file.to_xdmf(funcname, subnames, outputdir)
    file.close()
