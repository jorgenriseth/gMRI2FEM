from pathlib import Path
from typing import Callable, Optional

import dolfin as df
import nibabel
import numpy as np
import scipy
import simple_mri as sm
import pyvista as pv
from pantarei import FenicsStorage

from gmri2fem.utils import apply_affine, nan_filter_gaussian
from i2m.vtk2mri import mri_data_to_ndarray


def mri2fem_interpolate_quadrature(data_mri, V_target, quad_degree, mask=None):
    datamask = np.isfinite(data_mri.data)
    if mask is not None:
        datamask *= mask

    domain = V_target.mesh()
    quad_element = df.FiniteElement(
        "Quadrature", domain.ufl_cell(), degree=quad_degree, quad_scheme="default"
    )
    Q = df.FunctionSpace(domain, quad_element)
    dof_img_coordinates = locate_dof_voxels(Q, data_mri, rint=False)
    dof_img_voxels = find_dof_nearest_neighbours(dof_img_coordinates, datamask, N=1)

    dof_concentrations = data_mri.data[*dof_img_voxels]
    assert np.isfinite(dof_concentrations).all()
    q = df.Function(Q)
    q.vector()[:] = dof_concentrations

    dx = df.Measure("dx", metadata={"quadrature_degree": quad_degree})
    u, v = df.TrialFunction(V_target), df.TestFunction(V_target)
    a = u * v * dx
    l = q * v * dx
    A = df.assemble(a)
    b = df.assemble(l)
    uh = df.Function(V_target)
    df.solve(A, uh.vector(), b)
    return uh


def dolfin_mesh_to_pyvista_ugrid(mesh):
    # Extract the mesh coordinates and connectivity
    coordinates = mesh.coordinates()
    cells = mesh.cells()
    num_cells = cells.shape[0]
    num_points_per_cell = 4
    connectivity = np.hstack([np.full((num_cells, 1), num_points_per_cell), cells])
    return pv.UnstructuredGrid(
        connectivity, np.ones(num_cells) * pv.CellType.TETRA, coordinates
    )


def dolfin2mri(
    u: df.Function,
    reference_mri: sm.SimpleMRI,
    grid: pv.UnstructuredGrid = None,
    fieldname: str = "",
):
    ugrid = grid or dolfin_mesh_to_pyvista_ugrid(u.function_space().mesh())
    fname = fieldname or "u"
    ugrid.point_data[fname] = u.compute_vertex_values()

    shape = reference_mri.shape
    affine = reference_mri.affine
    image_grid = pv.ImageData(dimensions=shape).transform(affine, inplace=False)
    grid = image_grid.sample(ugrid, progress_bar=False, locator="static_cell")
    return sm.SimpleMRI(mri_data_to_ndarray(grid, fname, shape), affine)


def find_dof_nearest_neighbours(
    dof_inds: np.ndarray, mask: np.ndarray, N: int
) -> np.ndarray:
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


def mri2fem_interpolate_collocation(
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
    return mri2fem_interpolate_collocation(
        voxeldata, mri_volume.affine, functionspace, datafilter
    )


def fenicsstorage2xdmf(
    filepath, funcname: str, subnames: str | list[str], outputdir: Path
) -> None:
    file = FenicsStorage(filepath, "r")
    file.to_xdmf(funcname, subnames, outputdir)
    file.close()
