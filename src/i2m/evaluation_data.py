import click
import panta_rhei as pr
import simple_mri as sm
import dolfin as df
import numpy as np
import scipy
import tqdm
from petsc4py import PETSc
from dolfin import inner, grad
from glymphopt.utils import with_suffix


def scipy_sparse_csr_to_dolfin(A: scipy.sparse.csr_matrix) -> df.Matrix:
    petsc_mat = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr, A.indices, A.data))
    return df.Matrix(df.PETScMatrix(petsc_mat))


def create_evaluation_matrix(
    V: df.FunctionSpace,
    data_mris: list[sm.SimpleMRI],
):
    """
    ***NB:*** Although tempting to split this into one function for matrix creation,
    and one function for determination of valid voxel indices, this function
    needs to be one, since there is some memory bug occuring once some not
    yet known variable goes out of scope.
    """
    mesh = V.mesh()
    affine = data_mris[0].affine
    mask = np.prod([np.isfinite(mri.data) for mri in data_mris], axis=0)
    IJK = np.array(np.where(mask)).T
    XYZ = sm.apply_affine(affine, IJK)

    tree = mesh.bounding_box_tree()
    cells_containing_point = np.array(
        [tree.compute_first_entity_collision(df.Point(*xi)) for xi in tqdm.tqdm(XYZ)]
    )
    (points_of_interest,) = np.where(cells_containing_point < mesh.num_cells())
    cells_containing_point = cells_containing_point[points_of_interest]
    ijk = IJK[points_of_interest]
    xyz = XYZ[points_of_interest]
    M = scipy.sparse.lil_matrix((len(xyz), V.dim()))
    for i, (xi, cell_index_i) in enumerate(zip(tqdm.tqdm(xyz), cells_containing_point)):
        cell_global_dofs = V.dofmap().cell_dofs(cell_index_i)
        cell = df.Cell(V.mesh(), cell_index_i)
        dof_vals = V.element().evaluate_basis_all(
            xi, cell.get_vertex_coordinates(), cell.orientation()
        )
        M[i, cell_global_dofs] = dof_vals
    return M.tocsr(), ijk


@click.command("mesh_evaluation_data")
@click.argument("mris", type=click.Path(exists=True), nargs=-1)
@click.option("--input", "-i", type=click.Path(exists=True), required=True)
@click.option("--output", "-o", type=click.Path(), required=True)
@click.option("--csfmask", type=click.Path(exists=True))
@click.option("--timestamps", type=click.Path(exists=True))
def create_mesh_evaluation_data(mris, input, output, csfmask, timestamps: str):
    with df.HDF5File(df.MPI.comm_world, input, "r") as hdf:
        mesh = pr.read_domain(hdf)
    V = df.FunctionSpace(mesh, "CG", 1)
    data_mris = [sm.load_mri(p, dtype=np.single) for p in mris]
    if csfmask:
        mask_mri = sm.load_mri(csfmask, dtype=bool)
        for mri in data_mris:
            mri.data[mask_mri.data] = np.nan

    M, ijk = create_evaluation_matrix(V, data_mris)
    time = np.loadtxt(str(timestamps)) if timestamps else np.array([])
    np.savez_compressed(
        output,
        matrix_data=M.data,
        matrix_indices=M.indices,
        matrix_indptr=M.indptr,
        matrix_shape=M.shape,
        mri_indices=ijk,
        mri_shape=data_mris[0].shape,
        mri_affine=data_mris[0].affine.flatten(),
        **{f"vector{idx}": mri.data[*ijk.T] for idx, mri in enumerate(data_mris)},  # type: ignore
        time=time,
        fem_family="CG",
        fem_degree=1,
    )
