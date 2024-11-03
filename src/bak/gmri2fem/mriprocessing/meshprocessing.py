#!/usr/bin/env python
import logging
from pathlib import Path

import dolfin as df
import meshio
from dolfin import HDF5File, Mesh, MeshFunction, MeshValueCollection, XDMFFile
from pantarei import Domain

logger = logging.getLogger(__name__)


def mesh2xdmf(meshfile, xdmfdir):
    mesh = meshio.read(meshfile)
    logger.info(f"Converting {meshfile} to {xdmfdir}/xxxxx.xdmf")
    return meshio2xdmf(mesh, xdmfdir)


def meshio2xdmf(mesh, xdmfdir):
    polytope_label = "tetra"
    facet_label = "triangle"
    points = mesh.points
    polytopes = {polytope_label: mesh.cells_dict[polytope_label]}
    facets = {facet_label: mesh.cells_dict[facet_label]}

    meshdata = meshio.Mesh(points, polytopes)
    meshio.write("{}/mesh.xdmf".format(xdmfdir), meshdata)
    if "gmsh:physical" or "medit:ref" in mesh.cell_data_dict:
        cell_data_name = (
            "gmsh:physical" if "gmsh:physical" in mesh.cell_data_dict else "medit:ref"
        )
        # Write the subdomains of the mesh
        subdomains = {
            "subdomains": [mesh.cell_data_dict[cell_data_name][polytope_label]]
        }
        subdomainfile = meshio.Mesh(points, polytopes, cell_data=subdomains)
        meshio.write("{}/subdomains.xdmf".format(xdmfdir), subdomainfile)

        # Write the boundaries/interfaces of the mesh
        boundaries = {"boundaries": [mesh.cell_data_dict[cell_data_name][facet_label]]}
        boundaryfile = meshio.Mesh(points, facets, cell_data=boundaries)
        meshio.write("{}/boundaries.xdmf".format(xdmfdir), boundaryfile)

    return xdmfdir


def xdmf2hdf(xdmfdir, hdf5file):
    # Read xdmf-file into a FEniCS mesh
    dirpath = Path(xdmfdir)
    mesh = Mesh(df.MPI.comm_world)
    with XDMFFile(str(dirpath / "mesh.xdmf")) as meshfile:
        meshfile.read(mesh)

    # Read cell data to a MeshFunction (of dim n)
    n = mesh.topology().dim()
    subdomains = MeshFunction("size_t", mesh, n)
    with XDMFFile(str(dirpath / "subdomains.xdmf")) as subdomainfile:
        subdomainfile.read(subdomains, "subdomains")

    # Read face data into a Meshfunction of dim n-1
    bdrycollection = MeshValueCollection("size_t", mesh, n - 1)
    with XDMFFile(str(dirpath / "boundaries.xdmf")) as boundaryfile:
        boundaryfile.read(bdrycollection, "boundaries")
    boundaries = MeshFunction("size_t", mesh, bdrycollection)

    # Write all files into a single h5-file.
    with HDF5File(mesh.mpi_comm(), str(hdf5file), "w") as f:
        f.write(mesh, "/domain/mesh")
        f.write(subdomains, "/domain/subdomains")
        f.write(boundaries, "/domain/boundaries")
    return Path(hdf5file)


def hdf2fenics(hdf5file: Path, pack=False):
    """Function to read h5-file with annotated mesh, subdomains
    and boundaries into fenics mesh"""
    hdf = df.HDF5File(df.MPI.comm_world, str(hdf5file), "r")

    mesh = Mesh(df.MPI.comm_world)
    hdf.read(mesh, "/domain/mesh", False)
    n = mesh.topology().dim()
    subdomains = df.MeshFunction("size_t", mesh, n)
    boundaries = df.MeshFunction("size_t", mesh, n - 1)
    if hdf.has_dataset("/domain/subdomains"):
        hdf.read(subdomains, "/domain/subdomains")
    if hdf.has_dataset("/domain/boundaries"):
        hdf.read(boundaries, "/domain/boundaries")

    if pack:
        return Domain(mesh, subdomains, boundaries)

    return mesh, subdomains, boundaries
