from pathlib import Path

import click
import dolfin as df
import meshio
import numpy as np
import panta_rhei as pr
import pyvista as pv


def mesh_to_pyvista_grid(mesh: meshio.Mesh, dim: int):
    geometry = mesh_to_geometry(mesh, dim)
    cells = cell_matrix_to_vtk_cells(geometry["polytopes"]["tetra"])
    celltypes = [pv.CellType.TETRA] * len(geometry["polytopes"]["tetra"])
    points = geometry["points"]
    return pv.UnstructuredGrid(cells, celltypes, points)


def mesh_to_geometry(
    mesh: meshio.Mesh, dim: int
) -> dict[str, np.ndarray | dict[str, np.ndarray]]:
    if dim == 2:
        polytope_label = "triangle"
        facet_label = "line"
    elif dim == 3:
        polytope_label = "tetra"
        facet_label = "triangle"
    else:
        raise ValueError("dim should be in (2, 3), got {}.".format(dim))
    if dim == 2:
        points = mesh.points[:, :2]
    else:
        points = mesh.points
    polytopes = {polytope_label: mesh.cells_dict[polytope_label]}
    facets = {facet_label: mesh.cells_dict[facet_label]}
    return {"points": points, "polytopes": polytopes, "facets": facets}


def cell_matrix_to_vtk_cells(cells: np.ndarray):
    num_points_stack = cells.shape[1] * np.ones((cells.shape[0], 1), dtype=cells.dtype)
    return np.concatenate((num_points_stack, cells), axis=1).ravel()


def convert_to_vtk(hdf_data: Path, output: Path, ascii: bool = False):
    hdf = df.HDF5File(df.MPI.comm_world, str(hdf_data), "r")
    domain = pr.read_domain(hdf)
    concentration_times = pr.read_timevector(hdf, "concentration")
    boundary_times = pr.read_timevector(hdf, "boundary_concentration")
    concentrations = [
        pr.read_function(hdf, "concentration", domain, idx) for idx in range(5)
    ]
    boundary_concentrations = [
        pr.read_function(hdf, "boundary_concentration", domain, idx) for idx in range(5)
    ]
    dti = pr.read_function(hdf, "DTI", domain)
    md = pr.read_function(hdf, "MD", domain)
    fa = pr.read_function(hdf, "FA", domain)

    parcellations = df.MeshFunction("size_t", domain, 3)
    hdf.read(parcellations, "parcellations")
    pr.close(hdf)

    cells = cell_matrix_to_vtk_cells(domain.cells())
    celltypes = [pv.CellType.TETRA] * len(domain.cells())
    points = domain.coordinates()
    grid = pv.UnstructuredGrid(cells, celltypes, points)

    grid.cell_data["subdomains"] = domain.subdomains.array()
    grid.cell_data["parcellations"] = parcellations.array()
    for idx, (ci, cb_i) in enumerate(zip(concentrations, boundary_concentrations)):
        grid.point_data[f"concentration-{int(concentration_times[idx])}"] = (
            ci.compute_vertex_values()
        )
        grid.point_data[f"boundary-concentration-{int(boundary_times[idx])}"] = (
            cb_i.compute_vertex_values()
        )

    for name, func in [("dt", dti), ("fa", fa), ("md", md)]:
        output_vec = func.vector()[:]
        output_dim = output_vec.shape[0] // domain.cells().shape[0]
        grid.cell_data[name] = output_vec.reshape(-1, output_dim).squeeze()
    grid.save(output, binary=(not ascii))


@click.command()
@click.option("--input", "hdf_data", type=Path, required=True)
@click.option("--output", "output", type=Path, required=True)
@click.option("--ascii", "ascii", type=bool, is_flag=True)
def hdf2vtk(**kwargs):
    convert_to_vtk(**kwargs)


if __name__ == "__main__":
    hdf2vtk()
