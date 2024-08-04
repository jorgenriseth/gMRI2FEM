import json
import logging
from functools import partial
from pathlib import Path

import dolfin as df
import numpy as np
import pandas as pd
import pantarei as pr
import pint
from dolfin import grad, inner

from gmri2fem.models.parameters import fasttransfer_parameters, singlecomp_parameters
from gmri2fem.models.solvers import (
    ArtificialCSFConcentration,
    ArtificialSASConcentration,
    solve_diffusion,
)
from gmri2fem.utils import float_string_formatter

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
logging.getLogger("FFC").setLevel(logging.WARNING)
logging.getLogger("UFL").setLevel(logging.WARNING)
df.set_log_level(df.LogLevel.WARNING)
log = partial(pr.single_logger, logger)


def read_concentration_data(filepath, funcname) -> list[df.Function]:
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


def diffusion_model(
    coefficients: dict[str, float],
    inputfile: Path,
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

    # Setup timestepping
    dt = 3600  # s
    T = timevec[-1]
    time = pr.TimeKeeper(dt=dt, endtime=T)

    with df.HDF5File(df.MPI.comm_world, str(inputfile), "r") as f:
        if "source" in f.attributes("total_concentration"):
            f_code = f.attributes("total_concentration")["source"]
            f_coeffs = json.loads(f.attributes("total_concentration")["coeffs"])
            coefficients["source"] = df.Expression(f_code, **f_coeffs, t=time)

    if coefficients["robin"] == float("inf"):
        boundaries = [pr.DirichletBoundary(u_data, "everywhere")]
    elif coefficients["robin"] == 0:
        boundaries = []
    else:
        a = coefficients["robin"]
        log("info", "Using artificial data.")
        u_data = ArtificialCSFConcentration(
            sas=ArtificialSASConcentration(scale=0.52),
            ventricles=ArtificialSASConcentration(scale=0.2),
        )
        boundaries = [
            pr.RobinBoundary(a, u_data.outer, (4, 5)),
            pr.RobinBoundary(a, u_data.inner, 8),
        ]

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
    if args.model == "fasttransfer":
        defaults = fasttransfer_parameters(param_units)
    elif args.model == "singlecomp":
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
    parser.add_argument("--D", type=float)
    parser.add_argument("--r", type=float)
    parser.add_argument("--k", type=float)
    parser.add_argument("--visual", action="store_true")
    parser.add_argument("--noquant", action="store_true")
    parser.add_argument("--model", type=str, default="fasttransfer")
    args = parser.parse_args()

    defaults = set_default_coefficients(args.model)
    coefficients = {
        "D": float(args.D) if args.D is not None else defaults["D"],
        "r": float(args.r) if args.r is not None else defaults["r"],
        "robin": float(args.k) if args.k is not None else defaults["robin"],
    }
    if df.MPI.comm_world.rank == 0:
        from gmri2fem.models.parameters import print_quantities

        print()
        print("=== Coefficients: ===")
        print_quantities(coefficients, offset=10)
        print()
    computer = diffusion_model(coefficients, Path(args.input), Path(args.output))

    if args.visual:
        output = Path(args.output)
        file = pr.FenicsStorage(output, "r")
        k = float_string_formatter(float(coefficients["robin"]), 3)
        filename = f"visual/diffusion" + f"_robin" * (k != "inf") + ".xdmf"
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
