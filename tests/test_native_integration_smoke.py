"""Focused cross-package smoke tests for the Conda native dependency stack.

These are not a comprehensive feature suite. They cover representative
interactions that package resolution and import checks cannot validate:
FEniCS HDF5 I/O (fenics_io), SVMTK surface creation, and SimpleMRI-to-FEniCS
projection with both supported solver paths.
"""

import numpy as np
import dolfin as df
import SVMTK as svmtk


def test_fenics_io_hdf5_roundtrip_and_svmtk_surface(tmp_path):
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
    with df.HDF5File(mesh.mpi_comm(), path, "r") as stream:
        restored = read_domain(stream)
    assert restored.num_cells() == mesh.num_cells()
    surface = svmtk.Surface()
    surface.make_sphere(0, 0, 0, 1, 0.5)
    assert surface.num_vertices() > 0


def test_simple_mri_to_fenics_projection_preserves_constant_for_both_solvers():
    from simple_mri import SimpleMRI
    from i2m.mri2fenics import mri2fem_interpolate_quadrature

    V = df.FunctionSpace(df.UnitCubeMesh(2, 2, 2), "CG", 1)
    mri = SimpleMRI(np.ones((4, 4, 4), dtype=np.single), np.eye(4))
    for solver in ("direct", "iterative"):
        result = mri2fem_interpolate_quadrature(mri, V, 2, solver_type=solver)
        assert np.allclose(result.vector().get_local(), 1, atol=1e-6)
