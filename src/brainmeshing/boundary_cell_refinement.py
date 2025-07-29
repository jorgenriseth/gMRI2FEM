import click
import dolfin as df
import numpy as np
import pantarei as pr


@click.command("refine_mesh")
@click.option("--input", "-i", type=str, required=True)
@click.option("--output", "-o", type=str, required=True)
@click.option("--threshold", "-t", type=float, default=2.0)
@click.option("--function", "-f", "functions", type=str, multiple=True)
def refine_mesh_cli(input, output, threshold, functions):
    with df.HDF5File(df.MPI.comm_world, input, "r") as hdf:
        domain = pr.read_domain(hdf)
    refined_mesh = refine_mesh(domain, threshold)

    # TODO: Fix meshfunction reader.
    # TODO: Add facet-function mapping.
    for function in functions:
        subdomains = pr.read_function()
        map_mesh_function(subdomains, domain, refined_mesh)

    # TODO: Add writer to both HDF and XDMF.


def is_boundary_cell(cell: df.Cell) -> bool:
    return any(facet.exterior() for facet in df.facets(cell))


def refine_mesh(mesh, threshold):
    init_mesh = df.Mesh(mesh)
    potential_boundary_cells = df.MeshFunction(
        "size_t", init_mesh, init_mesh.topology().dim(), 1
    )

    while True:
        init_mesh.init()
        cells_marked_for_refinement = df.MeshFunction(
            "bool", init_mesh, init_mesh.topology().dim(), False
        )

        N = init_mesh.num_cells()
        is_large_boundary_cell = np.zeros(N, dtype=bool)
        num_cells_marked_for_refinement = 0
        for cell in filter(lambda c: potential_boundary_cells[c], df.cells(init_mesh)):  # type: ignore
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

    return init_mesh


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
