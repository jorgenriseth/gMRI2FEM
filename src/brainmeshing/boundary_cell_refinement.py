import click
import dolfin as df
import numpy as np
import pantarei as pr
import tqdm

from glymphopt.utils import with_suffix

df.parameters["refinement_algorithm"] = "plaza_with_parent_facets"


@click.command("refine_mesh")
@click.option("--input", "-i", type=str, required=True)
@click.option("--output", "-o", type=str, required=True)
@click.option("--threshold", "-t", type=float, default=2.0)
@click.option("--maxiter", type=int, default=5)
def refine_mesh_cli(input, output, threshold, maxiter):
    with df.HDF5File(df.MPI.comm_world, input, "r") as hdf:
        domain = pr.read_domain(hdf)
    refined_domain = refine_mesh(domain, threshold, max_iter=maxiter)

    with df.HDF5File(df.MPI.comm_world, output, "w") as hdf:
        pr.write_domain(hdf, refined_domain)

    with df.XDMFFile(
        df.MPI.comm_world, str(with_suffix(output, "_subdomains.xdmf"))
    ) as xdmf:
        xdmf.write(refined_domain.subdomains)

    with df.XDMFFile(
        df.MPI.comm_world, str(with_suffix(output, "_boundaries.xdmf"))
    ) as xdmf:
        xdmf.write(refined_domain.boundaries)


def is_boundary_cell(cell: df.Cell) -> bool:
    return any(facet.exterior() for facet in df.facets(cell))


def refine_mesh(domain, threshold, max_iter=5):
    init_mesh = df.Mesh(domain)
    potential_boundary_cells = df.MeshFunction(
        "size_t", init_mesh, init_mesh.topology().dim(), 1
    )
    init_subdomains = df.MeshFunction(
        "size_t", init_mesh, init_mesh.topology().dim(), 0
    )
    init_subdomains.array()[:] = domain.subdomains.array()
    init_boundaries = df.MeshFunction(
        "size_t", init_mesh, init_mesh.topology().dim() - 1, 0
    )
    init_boundaries.array()[:] = domain.boundaries.array()

    iter = 0
    while iter < max_iter:
        iter += 1
        init_mesh.init()
        cells_marked_for_refinement = df.MeshFunction(
            "bool", init_mesh, init_mesh.topology().dim(), False
        )
        N = init_mesh.num_cells()
        is_large_boundary_cell = np.zeros(N, dtype=bool)
        num_cells_marked_for_refinement = 0
        n = potential_boundary_cells.array()[:].sum()
        cell_generator = filter(
            lambda c: potential_boundary_cells[c],  # type: ignore
            df.cells(init_mesh),
        )
        for cell in tqdm.tqdm(cell_generator, total=n):
            if is_boundary_cell(cell) and cell.circumradius() > threshold:
                is_large_boundary_cell[cell.index()] = True

        potential_boundary_cells.array()[np.argwhere(~is_large_boundary_cell)] = 0  # type: ignore
        cells_marked_for_refinement.array()[:] = is_large_boundary_cell  # type: ignore
        num_cells_marked_for_refinement = is_large_boundary_cell.sum()
        print(f"Refinined {num_cells_marked_for_refinement} cells.")
        if num_cells_marked_for_refinement == 0:
            break
        refined_mesh = df.refine(init_mesh, cells_marked_for_refinement)
        potential_boundary_cells = df.adapt(potential_boundary_cells, refined_mesh)
        init_mesh = refined_mesh
        init_subdomains = df.adapt(init_subdomains, refined_mesh)
        init_boundaries = df.adapt(init_boundaries, refined_mesh)

    newfacets = init_boundaries.array() > len(init_boundaries.array())
    init_boundaries.array()[np.argwhere(newfacets)] = 0
    return pr.Domain(
        init_mesh,
        init_subdomains,
        init_boundaries,
    )


def map_mesh_function(source_mesh_func, source_mesh, dest_mesh):
    """
    Maps a cell-based MeshFunction from a source mesh to a destination mesh
    using DG-0 function interpolation.

    Args:
        source_mesh (df.Mesh): The mesh the source_mesh_func lives on.
        dest_mesh (df.Mesh): The mesh to map the function to.
        source_mesh_func (df.MeshFunction): The cell function to be mapped.

    Returns:
        df.MeshFunction: The mapped cell function on the destination mesh.
    """
    V_source = df.FunctionSpace(source_mesh, "DG", 0)
    dg_func_source = df.Function(V_source)
    dg_func_source.vector()[:] = source_mesh_func.array().astype(float)
    V_dest = df.FunctionSpace(dest_mesh, "DG", 0)
    dg_func_dest = df.Function(V_dest)
    dg_func_dest.interpolate(dg_func_source)
    dest_mesh_func = df.MeshFunction("size_t", dest_mesh, dest_mesh.topology().dim())
    dest_mesh_func.set_values(np.round(dg_func_dest.vector().get_local()).astype(int))
    return dest_mesh_func
