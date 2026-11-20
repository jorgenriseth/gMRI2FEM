"""Conversion of meshio-readable meshes (e.g. medit .mesh from SVMTK) to a
FEniCS HDF5 domain file, via an intermediate directory of XDMF files."""

from pathlib import Path

import dolfin as df
import meshio
import numpy as np

from fenics_io.domain import Domain
from fenics_io.io import write_domain

StrPath = str | Path

_TAG_KEYS = ("gmsh:physical", "medit:ref")


def mesh2xdmf(meshfile: StrPath, xdmfdir: StrPath, dim: int) -> Path:
    return meshio2xdmf(meshio.read(meshfile), xdmfdir, dim=dim)


def meshio2xdmf(mesh: meshio.Mesh, xdmfdir: StrPath, dim: int) -> Path:
    """Write `mesh.xdmf`, and if the mesh is tagged, `subdomains.xdmf` and
    `boundaries.xdmf` into `xdmfdir`."""
    if dim == 2:
        polytope_label, facet_label = "triangle", "line"
        points = mesh.points[:, :2]
    elif dim == 3:
        polytope_label, facet_label = "tetra", "triangle"
        points = mesh.points
    else:
        raise ValueError(f"dim should be in (2, 3), got {dim}.")

    xdmfdir = Path(xdmfdir)
    xdmfdir.mkdir(exist_ok=True, parents=True)
    polytopes = {polytope_label: mesh.cells_dict[polytope_label]}
    meshio.write(xdmfdir / "mesh.xdmf", meshio.Mesh(points, polytopes))

    tag_key = next((key for key in _TAG_KEYS if key in mesh.cell_data_dict), None)
    tags = mesh.cell_data_dict[tag_key] if tag_key is not None else {}
    written_tag_files = set()

    if polytope_label in tags:
        subdomains = {"subdomains": [tags[polytope_label]]}
        meshio.write(
            xdmfdir / "subdomains.xdmf",
            meshio.Mesh(points, polytopes, cell_data=subdomains),
        )
        written_tag_files.add("subdomains.xdmf")
    if facet_label in tags and facet_label in mesh.cells_dict:
        facets = {facet_label: mesh.cells_dict[facet_label]}
        boundaries = {"boundaries": [tags[facet_label]]}
        meshio.write(
            xdmfdir / "boundaries.xdmf",
            meshio.Mesh(points, facets, cell_data=boundaries),
        )
        written_tag_files.add("boundaries.xdmf")

    for filename in ("subdomains.xdmf", "boundaries.xdmf"):
        if filename not in written_tag_files:
            (xdmfdir / filename).unlink(missing_ok=True)
    return xdmfdir


def xdmf2hdf(xdmfdir: StrPath, hdf5file: StrPath) -> Path:
    """Read the XDMF files written by `meshio2xdmf` and store them as a Domain
    in `hdf5file`."""
    dirpath = Path(xdmfdir)
    mesh = df.Mesh()
    with df.XDMFFile(str(dirpath / "mesh.xdmf")) as meshfile:
        meshfile.read(mesh)

    n = mesh.topology().dim()
    subdomains = df.MeshFunction("size_t", mesh, n, 0)
    if (dirpath / "subdomains.xdmf").exists():
        with df.XDMFFile(str(dirpath / "subdomains.xdmf")) as subdomainfile:
            subdomainfile.read(subdomains, "subdomains")

    bdrycollection = df.MeshValueCollection("size_t", mesh, n - 1)
    if (dirpath / "boundaries.xdmf").exists():
        with df.XDMFFile(str(dirpath / "boundaries.xdmf")) as boundaryfile:
            boundaryfile.read(bdrycollection, "boundaries")
    boundaries = df.MeshFunction("size_t", mesh, bdrycollection)
    # Facets missing from the collection are set to max(size_t) by dolfin.
    untagged = boundaries.array() == np.iinfo(np.uintp).max
    boundaries.array()[untagged] = 0

    with df.HDF5File(mesh.mpi_comm(), str(hdf5file), "w") as f:
        write_domain(f, Domain(mesh, subdomains, boundaries))
    return Path(hdf5file)
