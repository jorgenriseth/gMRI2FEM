from typing import Optional
import click
import panta_rhei as pr
import simple_mri as sm
import dolfin as df
import numpy as np
import scipy
import tqdm
from dolfin import inner, grad
from glymphopt.utils import with_suffix


def petsc_to_scipy(petsc_mat):
    m, n = petsc_mat.getSize()
    indptr, indices, data = petsc_mat.getValuesCSR()
    return scipy.sparse.csr_matrix((data, indices, indptr), shape=(m, n))


def construct_minimal_evaluation_system(V, M):
    (nonzero_cols,) = np.where(M.getnnz(0) != 0)
    A = M.T @ M
    A = A[nonzero_cols].tocsr()

    u, v = df.TrialFunction(V), df.TestFunction(V)
    dx = df.Measure("dx", domain=V.mesh())
    K = df.assemble(inner(grad(u), grad(v)) * dx)

    K_sp = petsc_to_scipy(df.as_backend_type(K).mat())
    AA = scipy.sparse.bmat(
        [[K_sp, A.T], [A, scipy.sparse.csr_matrix((A.shape[0], A.shape[0]))]]
    ).tocsr()
    return AA


def construct_minimal_evaluation_vector(M, x):
    (nonzero_cols,) = np.where(M.getnnz(0) != 0)
    b = M.T @ x
    bb = np.zeros(len(b) + len(nonzero_cols))
    bb[len(b) :] = b[nonzero_cols]
    return bb


def map_mri_to_mesh(
    C: list[np.ndarray], V: df.FunctionSpace, M: scipy.sparse.csr_matrix, verbose=True
) -> list[df.Function]:
    A = construct_minimal_evaluation_system(V, M)
    B = [construct_minimal_evaluation_vector(M, c) for c in C]
    U = [df.Function(V, name="concentration") for _ in range(len(C))]
    for idx, b in enumerate(tqdm.tqdm(B, disable=(not verbose))):
        uu, returncode = scipy.sparse.linalg.minres(A, b, show=verbose)
        U[idx].vector()[:] = uu[: V.dim()]
    return U


def reconstruct_evaluation_matrix(
    npzfile: np.lib.npyio.NpzFile,
) -> scipy.sparse.csr_matrix:
    return scipy.sparse.csr_matrix(
        (
            npzfile["matrix_data"],
            npzfile["matrix_indices"],
            npzfile["matrix_indptr"],
        ),
        shape=npzfile["matrix_shape"],
    )


@click.command()
@click.option("--mesh", "-m", type=str, required=True)
@click.option("--eval", "-e", type=str, required=True)
@click.option("--output", "-o", type=str, required=True)
@click.option("--timestamps", type=click.Path(exists=True))
@click.option("--verbose", type=bool, is_flag=True)
@click.option("--visual", type=bool, is_flag=True)
def map_evaluation_data_to_mesh(
    mesh: str,
    eval: str,
    output: str,
    timestamps: Optional[str],
    verbose: bool = True,
    visual: bool = True,
):
    with df.HDF5File(df.MPI.comm_world, mesh, "r") as hdf:
        mesh = pr.read_domain(hdf)

    npzfile = np.load(eval)
    M = reconstruct_evaluation_matrix(npzfile)
    # scipy.sparse.csr_matrix(
    #     (
    #         npzfile["matrix_data"],
    #         npzfile["matrix_indices"],
    #         npzfile["matrix_indptr"],
    #     ),
    #     shape=npzfile["matrix_shape"],
    # )
    C = [npzfile[label] for idx in range(5) if (label := f"vector{idx}") in npzfile]
    assert len(C) > 0, f"No vector data found in {eval}"

    if timestamps:
        time = np.loadtxt(timestamps)
        assert len(time) == len(C), (
            f"Length of timestamps in {timestamps} does not match length of data"
        )
    else:
        time = npzfile["time"]
        if len(time) == 0:
            raise RuntimeError(f"Timestamps not provided, and not found in {eval}")

    family, degree = str(npzfile["fem_family"]), int(npzfile["fem_degree"])

    V = df.FunctionSpace(mesh, family, degree)
    U = map_mri_to_mesh(C, V, M, verbose)

    with df.HDF5File(df.MPI.comm_world, output, "w") as hdf:
        pr.write_domain(hdf, mesh)
        for t, u in zip(time, U):
            pr.write_checkpoint(hdf, u, str(u.name()), t=t)

    if visual:
        with df.XDMFFile(df.MPI.comm_world, str(with_suffix(output, ".xdmf"))) as xdmf:
            for t, u in zip(time, U):
                xdmf.write(u, t=t)
