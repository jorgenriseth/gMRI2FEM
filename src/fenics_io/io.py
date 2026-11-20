"""Reading and writing FEniCS meshes and (time-dependent) functions to HDF5.

File layout:
    /domain/mesh, /domain/subdomains, /domain/boundaries
    /<name>/vector_<idx>  (with attribute "timestamp"), and attributes
    "signature" and "count" on /<name>, as written by dolfin.HDF5File.
"""

import ast
from pathlib import Path
from typing import Callable, Optional

import dolfin as df
import numpy as np
import ufl_legacy
from loguru import logger

from fenics_io.domain import Domain

StrPath = str | Path


def open_hdf(filepath: StrPath, mode: str) -> df.HDF5File:
    """Open an HDF5File on MPI.comm_world. In write/append mode the parent
    directory is created, and append falls back to write if the file is missing."""
    path = Path(filepath)
    comm = df.MPI.comm_world
    if mode in ("w", "a"):
        if comm.rank == 0:
            path.parent.mkdir(exist_ok=True, parents=True)
        comm.barrier()
        if mode == "a" and not path.exists():
            mode = "w"
    return df.HDF5File(comm, str(path), mode)


def write_domain(hdf: df.HDF5File, domain: df.Mesh) -> None:
    hdf.write(domain, "/domain/mesh")
    if isinstance(domain, Domain):
        hdf.write(domain.subdomains, "/domain/subdomains")
        hdf.write(domain.boundaries, "/domain/boundaries")


def read_domain(hdf: df.HDF5File | StrPath) -> Domain:
    """Read a Domain from an open HDF5File or from a file path."""
    if not isinstance(hdf, df.HDF5File):
        with df.HDF5File(df.MPI.comm_world, str(hdf), "r") as f:
            return read_domain(f)

    mesh = df.Mesh(df.MPI.comm_world)
    hdf.read(mesh, "/domain/mesh", False)
    n = mesh.topology().dim()
    subdomains = df.MeshFunction("size_t", mesh, n, 0)
    boundaries = df.MeshFunction("size_t", mesh, n - 1, 0)
    if hdf.has_dataset("/domain/subdomains"):
        hdf.read(subdomains, "/domain/subdomains")
    if hdf.has_dataset("/domain/boundaries"):
        hdf.read(boundaries, "/domain/boundaries")
    return Domain(mesh, subdomains, boundaries)


def read_element(hdf: df.HDF5File, function_name: str) -> ufl_legacy.FiniteElementBase:
    signature = hdf.attributes(function_name)["signature"]
    return _read_signature(signature)


def write_function(hdf: df.HDF5File, function: df.Function, name: str) -> None:
    """Write the first entry of a function series (t=0.0), together with the
    mesh if no domain is stored in the file yet."""
    if not hdf.has_dataset("/domain"):
        write_domain(hdf, function.function_space().mesh())
    hdf.write(function, name, 0.0)


def read_function(
    hdf: df.HDF5File, name: str, domain: Optional[df.Mesh] = None, idx: int = 0
) -> df.Function:
    if domain is None:
        domain = read_domain(hdf)
    element = read_element(hdf, name)
    V = df.FunctionSpace(domain, element)
    u = df.Function(V, name=name)
    hdf.read(u, f"{name}/vector_{idx}")
    return u


def write_checkpoint(
    hdf: df.HDF5File, function: df.Function, name: str, t: float
) -> None:
    hdf.write(function, name, t)


def read_checkpoint(
    hdf: df.HDF5File, u: df.Function, name: str, idx: int
) -> df.Function:
    hdf.read(u, f"{name}/vector_{idx}")
    return u


def read_checkpoint_time(hdf: df.HDF5File, name: str, idx: int) -> float:
    return hdf.attributes(f"{name}/vector_{idx}")["timestamp"]


def read_timevector(hdf: df.HDF5File, function_name: str) -> np.ndarray:
    num_entries = int(hdf.attributes(function_name)["count"])
    return np.array(
        [read_checkpoint_time(hdf, function_name, i) for i in range(num_entries)]
    )


def hdf_to_xdmf(
    hdf: df.HDF5File,
    funcname: str,
    subnames: str | list,
    outputpattern: StrPath | Callable[[str], Path],
) -> None:
    """Write all checkpoints of `funcname` to XDMF. For mixed functions,
    `subnames` is a (nested) list naming each subfunction, each written to the
    file `outputpattern(subname)`. If `outputpattern` is a directory, files are
    named `<dir>/<funcname>_<subname>.xdmf`."""
    if not callable(outputpattern):
        outdir = Path(outputpattern)
        outputpattern = lambda x: outdir / f"{funcname}_{x}.xdmf"

    xdmfs = {}
    try:
        for name in _flat(subnames):
            filename = Path(outputpattern(name))
            filename.parent.mkdir(exist_ok=True, parents=True)
            xdmfs[name] = df.XDMFFile(df.MPI.comm_world, str(filename))

        times = read_timevector(hdf, funcname)
        u = read_function(hdf, funcname)
        for idx, ti in enumerate(times):
            logger.debug(f"Writing {funcname} checkpoint {idx + 1}/{len(times)}")
            u = read_checkpoint(hdf, u, funcname, idx)
            _write_to_xdmf(xdmfs, u, ti, subnames)
    finally:
        for xdmf in xdmfs.values():
            xdmf.close()


def _write_to_xdmf(xdmfs, u, t, names):
    if isinstance(names, str):
        u.rename(names, "")
        xdmfs[names].write(u, t)
    else:
        for uj, name in zip(u.split(deepcopy=True), names):
            _write_to_xdmf(xdmfs, uj, t, name)


def _flat(pool) -> list[str]:
    if isinstance(pool, str):
        return [pool]
    res = []
    for v in pool:
        res += _flat(v)
    return res


def _read_signature(signature: str) -> ufl_legacy.FiniteElementBase:
    """Reconstruct a UFL element from DOLFIN's serialized element signature."""
    from ufl_legacy import Cell

    elements = {
        "FiniteElement": df.FiniteElement,
        "MixedElement": df.MixedElement,
        "TensorElement": df.TensorElement,
        "VectorElement": df.VectorElement,
    }
    constructors = {
        **elements,
        "Cell": Cell,
    }
    cells = {
        "interval": df.interval,
        "triangle": df.triangle,
        "tetrahedron": df.tetrahedron,
        "quadrilateral": df.quadrilateral,
        "hexahedron": df.hexahedron,
    }

    try:
        expression = ast.parse(signature, mode="eval").body
    except SyntaxError as exc:
        raise ValueError("Invalid UFL element signature") from exc

    if (
        not isinstance(expression, ast.Call)
        or not isinstance(expression.func, ast.Name)
        or expression.func.id not in elements
    ):
        raise ValueError("Signature must call a supported UFL element constructor")

    def parse_value(node: ast.expr):
        if isinstance(node, ast.Constant):
            if node.value is None or type(node.value) in (
                str,
                int,
                float,
                bool,
            ):
                return node.value
        elif isinstance(node, ast.Name) and node.id in cells:
            return cells[node.id]
        elif isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name) or node.func.id not in constructors:
                raise ValueError("Unsupported constructor in UFL element signature")
            if any(isinstance(arg, ast.Starred) for arg in node.args):
                raise ValueError("Starred arguments are not allowed in UFL signatures")
            if any(keyword.arg is None for keyword in node.keywords):
                raise ValueError("Expanded keyword arguments are not allowed in UFL signatures")
            args = [parse_value(arg) for arg in node.args]
            kwargs = {
                keyword.arg: parse_value(keyword.value) for keyword in node.keywords
            }
            return constructors[node.func.id](*args, **kwargs)
        elif isinstance(node, (ast.List, ast.Tuple, ast.Set)):
            values = [parse_value(item) for item in node.elts]
            if isinstance(node, ast.List):
                return values
            if isinstance(node, ast.Tuple):
                return tuple(values)
            return set(values)
        elif isinstance(node, ast.Dict):
            if any(key is None for key in node.keys):
                raise ValueError("Dictionary unpacking is not allowed in UFL signatures")
            return {
                parse_value(key): parse_value(value)
                for key, value in zip(node.keys, node.values)
            }
        elif isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = parse_value(node.operand)
            if type(value) in (int, float):
                return value if isinstance(node.op, ast.UAdd) else -value
        raise ValueError("Unsupported expression in UFL element signature")

    return parse_value(expression)
