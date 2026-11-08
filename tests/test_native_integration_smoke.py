"""Focused cross-package smoke tests for the Conda native dependency stack.

These are not a comprehensive feature suite. They cover representative
interactions that package resolution and import checks cannot validate:
dolfin HDF5 I/O through fenics_io, SVMTK surface creation, meshio XDMF
(h5py) to dolfin mesh conversion, and SimpleMRI-to-FEniCS projection with
both supported solver paths.
"""

import numpy as np
import dolfin as df


def test_fenics_io_hdf5_roundtrip(tmp_path):
    from fenics_io import Domain, read_domain, write_domain

    mesh = df.UnitCubeMesh(2, 2, 2)
    domain = Domain(
        mesh,
        df.MeshFunction("size_t", mesh, 3, 1),
        df.MeshFunction("size_t", mesh, 2, 0),
    )
    path = str(tmp_path / "mesh.hdf")
    with df.HDF5File(mesh.mpi_comm(), path, "w") as stream:
        write_domain(stream, domain)
    restored = read_domain(path)
    assert restored.num_cells() == mesh.num_cells()
    assert (restored.subdomains.array() == 1).all()


def test_svmtk_surface():
    import SVMTK as svmtk

    surface = svmtk.Surface()
    surface.make_sphere(0, 0, 0, 1, 0.5)
    assert surface.num_vertices() > 0


def test_meshio_xdmf_to_fenics_domain(tmp_path):
    import meshio
    from fenics_io import meshio2xdmf, read_domain, xdmf2hdf

    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    mesh = meshio.Mesh(
        points,
        [("tetra", np.array([[0, 1, 2, 3]])), ("triangle", np.array([[0, 1, 2]]))],
        cell_data={"medit:ref": [np.array([3]), np.array([4])]},
    )
    meshio2xdmf(mesh, tmp_path / "xdmf", dim=3)
    domain = read_domain(xdmf2hdf(tmp_path / "xdmf", tmp_path / "mesh.hdf"))
    assert list(domain.subdomains.array()) == [3]
    assert 4 in domain.boundaries.array()


def test_simple_mri_to_fenics_projection_preserves_constant_for_both_solvers():
    from simple_mri import SimpleMRI
    from i2m.mri2fenics import mri2fem_interpolate_quadrature

    V = df.FunctionSpace(df.UnitCubeMesh(2, 2, 2), "CG", 1)
    mri = SimpleMRI(np.ones((4, 4, 4), dtype=np.single), np.eye(4))
    for solver in ("direct", "iterative"):
        result = mri2fem_interpolate_quadrature(mri, V, 2, solver_type=solver)
        assert np.allclose(result.vector().get_local(), 1, atol=1e-6)
