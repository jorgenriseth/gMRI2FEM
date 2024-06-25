from functools import partial
from pathlib import Path
from typing import Optional

import dolfin as df
import numpy as np
import pandas as pd
import pantarei as pr
from dolfin import grad, inner
from loguru import logger

from gmri2fem.models.parameters import fasttransfer_parameters, singlecomp_parameters
from gmri2fem.models.solvers import solve_diffusion
from gmri2fem.utils import float_string_formatter


def read_concentration_data(filepath, funcname) -> tuple[np.ndarray, list[df.Function]]:
    store = pr.FenicsStorage(filepath, "r")
    tvec = store.read_timevector(funcname)
    c = store.read_function(funcname, idx=0)
    C = [df.Function(c.function_space()) for _ in range(tvec.size)]
    for idx in range(len(C)):
        store.read_checkpoint(C[idx], funcname, idx)
    return tvec, C


def diffusion_form(V, coefficients, boundaries, u0, dt):
    u, v = df.TrialFunction(V), df.TestFunction(V)
    D = coefficients["D"]
    r = coefficients["r"]
    f = coefficients["source"] if "source" in coefficients else 0
    dx = df.Measure("dx", domain=V.mesh())
    F = (((u - u0) / dt + (r * u - f)) * v + inner(D * grad(u), grad(v))) * dx
    return F + pr.process_boundary_forms(u, v, boundaries)


def read_dti(dti_path: Path, mesh: df.Mesh, funcname: str = "DTI") -> df.Function:
    hdf = df.HDF5File(df.MPI.comm_world, str(dti_path), "r")
    D = pr.read_function(hdf, funcname, domain=mesh)
    hdf.close()
    return D


def diffusion_model_with_dti(
    coefficients: dict[str, float],
    inputfile: Path,
    dti_path: Path,
    output: Path,
    method="gmres",
    preconditioner="hypre_amg",
):
    # Read concentration data for Dirichlet BC.
    timevec, data = read_concentration_data(inputfile, "total_concentration")
    u_data = pr.DataInterpolator(data, timevec)
    u0 = u_data.copy(deepcopy=True)

    # Read functionspace and domain from
    V = u_data.function_space()
    domain = V.mesh()
    coefficients["D"] = read_dti(dti_path, domain)

    # Setup timestepping
    dt = 3600  # s
    T = timevec[-1]
    time = pr.TimeKeeper(dt=dt, endtime=T)

    boundaries = [pr.DirichletBoundary(u_data, "everywhere")]
    dx = df.Measure("dx", domain=domain, subdomain_data=domain.subdomains)
    computer = solve_diffusion(
        u_data=u_data,
        u0=u0,
        V=V,
        form=diffusion_form,
        coefficients=coefficients,
        boundaries=boundaries,
        time=time,
        solver=pr.StationaryProblemSolver(method, preconditioner),
        storage=pr.FenicsStorage(output, "w"),
        computer=solute_quantifier(dx),
    )
    return computer


def solute_quantifier(dx):
    return pr.BaseComputer(
        {
            "whole-brain": lambda u: df.assemble(u * dx),
            "gray-matter": lambda u: df.assemble(u * dx(1)),
            "white-matter": lambda u: df.assemble(u * dx(2)),
        }
    )


def set_default_coefficients(model: str):
    param_units = {"D": "mm**2 / s", "r": "1 / s", "robin": "mm / s"}
    if model == "fasttransfer":
        defaults = fasttransfer_parameters(param_units)
    elif model == "singlecomp":
        defaults = singlecomp_parameters(param_units)
    else:
        raise ValueError(
            f"Invalid model '{model}' should be 'fasttransfer' or 'singlecomp'."
        )
    return defaults


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input data file", type=Path, required=True)
    parser.add_argument(
        "--output", type=Path, help="Filepath of output .hdf", required=True
    )
    parser.add_argument("--dti_data", type=Path, required=True)
    parser.add_argument("--visual", action="store_true")
    parser.add_argument("--noquant", action="store_true")
    args = parser.parse_args()

    defaults = set_default_coefficients("singlecomp")
    coefficients = {"D": float(0.0), "r": float(defaults["r"]), "robin": float("inf")}
    if df.MPI.comm_world.rank == 0:
        from gmri2fem.models.parameters import print_quantities

        print()
        print("=== Coefficients: ===")
        print_quantities(coefficients, offset=10)
        print()

    computer = diffusion_model_with_dti(
        coefficients, args.input, args.dti_data, args.output
    )

    if args.visual:
        output = Path(args.output)
        file = pr.FenicsStorage(output, "r")
        k = float_string_formatter(float(coefficients["robin"]), 3)
        filename = f"visual/diffusion" + "_robin" * (k != "inf") + ".xdmf"
        file.to_xdmf(
            "total_concentration",
            "total_concentration",
            lambda _: output.parent / filename,
        )
        file.close()

    if df.MPI.comm_world.rank == 0 and not args.noquant:
        logger.info("Building dataframe from computer.")
        dframe = pd.DataFrame.from_dict(computer.values)
        logger.info("Storing dataframe to csv")
        dframe.to_csv(Path(args.output).with_suffix(".csv"))
