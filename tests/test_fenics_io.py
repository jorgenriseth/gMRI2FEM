import dolfin as df
import meshio
import numpy as np
import pytest

import fenics_io as fio


def _tagged_domain():
    mesh = df.UnitCubeMesh(2, 2, 2)
    subdomains = df.MeshFunction("size_t", mesh, 3, 0)
    subdomains.array()[:] = np.arange(mesh.num_cells()) % 3 + 1
    boundaries = df.MeshFunction("size_t", mesh, 2, 0)
    df.CompiledSubDomain("on_boundary && near(x[0], 0)").mark(boundaries, 7)
    return fio.Domain(mesh, subdomains, boundaries)


def _assert_same_domain(a, b):
    assert a.num_cells() == b.num_cells()
    assert np.allclose(a.coordinates(), b.coordinates())
    assert np.array_equal(a.subdomains.array(), b.subdomains.array())
    assert np.array_equal(a.boundaries.array(), b.boundaries.array())


def test_domain_roundtrip_from_handle_and_path(tmp_path):
    domain = _tagged_domain()
    path = tmp_path / "sub" / "domain.hdf"
    with fio.open_hdf(path, "w") as hdf:
        fio.write_domain(hdf, domain)
    with df.HDF5File(df.MPI.comm_world, str(path), "r") as hdf:
        _assert_same_domain(domain, fio.read_domain(hdf))
    _assert_same_domain(domain, fio.read_domain(path))
    _assert_same_domain(domain, fio.read_domain(str(path)))


def test_read_domain_without_tags_defaults_to_zero(tmp_path):
    mesh = df.UnitSquareMesh(2, 2)
    path = tmp_path / "mesh.hdf"
    with fio.open_hdf(path, "w") as hdf:
        fio.write_domain(hdf, mesh)
    domain = fio.read_domain(path)
    assert (domain.subdomains.array() == 0).all()
    assert (domain.boundaries.array() == 0).all()


def test_open_hdf_append_creates_missing_file(tmp_path):
    path = tmp_path / "new" / "file.hdf"
    with fio.open_hdf(path, "a") as hdf:
        fio.write_domain(hdf, df.UnitIntervalMesh(3))
    assert path.exists()


def test_read_signature_supports_cell_constructor():
    element = fio.io._read_signature("FiniteElement('CG', Cell('triangle'), 1)")
    assert element == df.FiniteElement("CG", df.triangle, 1)


@pytest.mark.parametrize(
    "make_space",
    [
        lambda m: df.FunctionSpace(m, "CG", 1),
        lambda m: df.VectorFunctionSpace(m, "CG", 2),
        lambda m: df.TensorFunctionSpace(m, "DG", 0),
        lambda m: df.FunctionSpace(
            m,
            df.MixedElement(
                [df.FiniteElement("CG", m.ufl_cell(), 1)] * 2
            ),
        ),
    ],
    ids=["CG1", "vectorCG2", "tensorDG0", "mixed"],
)
def test_function_checkpoint_roundtrip(tmp_path, make_space):
    domain = _tagged_domain()
    V = make_space(domain)
    times = [0.0, 1.5, 3.0]
    values = []
    path = tmp_path / "data.hdf"
    with fio.open_hdf(path, "w") as hdf:
        for idx, t in enumerate(times):
            u = df.Function(V)
            u.vector()[:] = np.random.default_rng(idx).random(V.dim())
            values.append(u.vector().get_local())
            if idx == 0:
                fio.write_function(hdf, u, "u")
            else:
                fio.write_checkpoint(hdf, u, "u", t)

    with df.HDF5File(df.MPI.comm_world, str(path), "r") as hdf:
        assert np.allclose(fio.read_timevector(hdf, "u"), times)
        assert fio.read_checkpoint_time(hdf, "u", 1) == pytest.approx(1.5)
        assert fio.read_element(hdf, "u") == V.ufl_element()
        u = fio.read_function(hdf, "u")
        assert np.allclose(u.vector().get_local(), values[0])
        _assert_same_domain(domain, u.function_space().mesh())
        for idx in range(len(times)):
            fio.read_checkpoint(hdf, u, "u", idx)
            assert np.allclose(u.vector().get_local(), values[idx])
        u2 = fio.read_function(hdf, "u", domain, idx=2)
        assert np.allclose(u2.vector().get_local(), values[2])


@pytest.mark.parametrize(
    "signature",
    [
        "__import__('os')",
        "FiniteElement.__class__",
        "FiniteElement('CG', triangle, [x for x in ()])",
        "FiniteElement('CG', triangle, 1 + 1)",
    ],
)
def test_read_signature_rejects_executable_expressions(signature):
    with pytest.raises(ValueError):
        fio.io._read_signature(signature)


def test_hdf_to_xdmf(tmp_path):
    mesh = df.UnitSquareMesh(2, 2)
    V = df.FunctionSpace(mesh, "CG", 1)
    W = df.FunctionSpace(
        mesh, df.MixedElement([df.FiniteElement("CG", mesh.ufl_cell(), 1)] * 2)
    )
    path = tmp_path / "data.hdf"
    with fio.open_hdf(path, "w") as hdf:
        fio.write_function(hdf, df.Function(V), "c")
        fio.write_checkpoint(hdf, df.Function(V), "c", 1.0)
        fio.write_function(hdf, df.Function(W), "w")

    with df.HDF5File(df.MPI.comm_world, str(path), "r") as hdf:
        fio.hdf_to_xdmf(hdf, "c", "internal", lambda _: tmp_path / "out" / "c.xdmf")
        fio.hdf_to_xdmf(hdf, "w", ["a", "b"], tmp_path / "mixed")
    assert (tmp_path / "out" / "c.xdmf").exists()
    assert (tmp_path / "mixed" / "w_a.xdmf").exists()
    assert (tmp_path / "mixed" / "w_b.xdmf").exists()


def _medit_mesh(tagged: bool) -> meshio.Mesh:
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float
    )
    tetra = np.array([[0, 1, 2, 3], [1, 2, 3, 4]])
    tris = np.array([[0, 1, 2], [1, 2, 4]])
    cells = [("tetra", tetra), ("triangle", tris)]
    if not tagged:
        return meshio.Mesh(points, cells)
    return meshio.Mesh(
        points, cells, cell_data={"medit:ref": [np.array([1, 2]), np.array([5, 6])]}
    )


def test_mesh2xdmf_xdmf2hdf(tmp_path):
    meshfile = tmp_path / "brain.mesh"
    meshio.write(meshfile, _medit_mesh(tagged=True))
    xdmfdir = tmp_path / "xdmfs"
    fio.mesh2xdmf(meshfile, xdmfdir, dim=3)
    assert (xdmfdir / "subdomains.xdmf").exists()
    assert (xdmfdir / "boundaries.xdmf").exists()

    domain = fio.read_domain(fio.xdmf2hdf(xdmfdir, tmp_path / "brain.hdf"))
    assert domain.num_cells() == 2
    assert sorted(domain.subdomains.array()) == [1, 2]
    assert sorted(set(domain.boundaries.array())) == [0, 5, 6]


def test_meshio2xdmf_without_tags(tmp_path):
    xdmfdir = tmp_path / "xdmfs"
    fio.meshio2xdmf(_medit_mesh(tagged=True), xdmfdir, dim=3)
    assert (xdmfdir / "subdomains.xdmf").exists()
    assert (xdmfdir / "boundaries.xdmf").exists()

    fio.meshio2xdmf(_medit_mesh(tagged=False), xdmfdir, dim=3)
    assert (xdmfdir / "mesh.xdmf").exists()
    assert not (xdmfdir / "subdomains.xdmf").exists()
    assert not (xdmfdir / "boundaries.xdmf").exists()

    domain = fio.read_domain(fio.xdmf2hdf(xdmfdir, tmp_path / "brain.hdf"))
    assert domain.num_cells() == 2
    assert (domain.subdomains.array() == 0).all()
    assert (domain.boundaries.array() == 0).all()


def test_compatible_with_panta_rhei_files(tmp_path):
    pr = pytest.importorskip("panta_rhei")
    domain = _tagged_domain()
    V = df.FunctionSpace(domain, "CG", 1)
    u = df.interpolate(df.Expression("x[0]", degree=1), V)

    legacy = tmp_path / "legacy.hdf"
    with df.HDF5File(df.MPI.comm_world, str(legacy), "w") as hdf:
        pr.write_domain(hdf, pr.Domain(domain, domain.subdomains, domain.boundaries))
        pr.write_function(hdf, u, "u")
        pr.write_checkpoint(hdf, u, "u", 2.0)
    with df.HDF5File(df.MPI.comm_world, str(legacy), "r") as hdf:
        _assert_same_domain(domain, fio.read_domain(hdf))
        assert np.allclose(fio.read_timevector(hdf, "u"), [0.0, 2.0])
        v = fio.read_function(hdf, "u", idx=1)
        assert np.allclose(v.vector().get_local(), u.vector().get_local())

    new = tmp_path / "new.hdf"
    with fio.open_hdf(new, "w") as hdf:
        fio.write_domain(hdf, domain)
        fio.write_function(hdf, u, "u")
    with df.HDF5File(df.MPI.comm_world, str(new), "r") as hdf:
        legacy_domain = pr.read_domain(hdf)
        assert np.array_equal(
            legacy_domain.subdomains.array(), domain.subdomains.array()
        )
        v = pr.read_function(hdf, "u")
        assert np.allclose(v.vector().get_local(), u.vector().get_local())
