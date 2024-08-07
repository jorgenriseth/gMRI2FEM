import logging
import time as pytime
from functools import partial
from typing import Callable, Optional

import dolfin as df
import pantarei as pr
import numpy as np

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
logging.getLogger("FFC").setLevel(logging.WARNING)
logging.getLogger("UFL").setLevel(logging.WARNING)
df.set_log_level(df.LogLevel.WARNING)
log = partial(pr.single_logger, logger)


class ArtificialSASConcentration(df.Constant):
    def __init__(self, scale=0.25, t1=4.43e4, t2=8.5e4):
        super().__init__(0)
        self.s = scale
        self.t1 = t1
        self.t2 = t2

    def update(self, t: float) -> df.Function:
        newval = self.s * (-np.exp(-t / self.t1) + np.exp(-t / self.t2))
        self.assign(newval)

class ArtificialCSFConcentration:
    def __init__(self, sas, ventricles):
        self.inner = ventricles
        self.outer = sas
            
    def update(self, t: float) -> df.Function:
        self.inner.update(t)
        self.outer.update(t)

def repeated_assigner(u: df.Function, ui: df.Function):
    """Assigns to each component of a function u - residing in a vector function
    space W - a function ui residing in the component subspace V."""
    W = u.function_space()
    V = ui.function_space()
    n = W.num_sub_spaces()
    df.FunctionAssigner(W, [V] * n).assign(u, [ui] * n)
    return u


def solve_multidiffusion(
    u_data: pr.DataInterpolator,
    W: df.FunctionSpace,
    form: pr.TimedependentForm,
    coefficients: pr.CoefficientsDict,
    boundaries: list[pr.BoundaryData],
    time: pr.TimeKeeper,
    solver: pr.StationaryProblemSolver,
    storage: pr.FenicsStorage,
    total_storage: pr.FenicsStorage,
    total: Callable[[df.Function], df.Function],
    computer: Optional[pr.BaseComputer] = None,
):
    computer = pr.set_optional(computer, pr.NullComputer)

    # Define boundayr condition, all zeros if artificial data.
    if isinstance(u_data, (ArtificialSASConcentration, ArtificialCSFConcentration)):
        u0 = df.Function(W)
    else: 
        u0 = repeated_assigner(df.Function(W), u_data)

    storage.write_function(u0, "fluid_concentrations")
    uT = total(u0)
    computer.compute(time, uT)
    total_storage.write_function(uT, "total_concentration")

    dirichlet_bcs = pr.process_dirichlet(W, boundaries)

    F = form(W, coefficients, boundaries, u0, time.dt)
    a = df.lhs(F)
    l = df.rhs(F)
    A = df.assemble(a)

    u = df.Function(W, name="fluid_concentrations")
    tic = pytime.time()
    for ti in time:
        pr.print_progress(float(ti), time.endtime, rank=df.MPI.comm_world.rank)
        u_data.update(float(ti))
        b = df.assemble(l)
        solver.solve(u, A, b, dirichlet_bcs)

        uT = total(u)
        computer.compute(ti, uT)
        storage.write_checkpoint(u, "fluid_concentrations", float(ti))
        total_storage.write_checkpoint(uT, "total_concentration", float(ti))
        u0.assign(u)

    storage.close()
    total_storage.close()

    log("info", "Time loop finished.")
    toc = pytime.time()
    df.MPI.comm_world.barrier()
    log("info", f"Elapsed time in loop: {toc - tic:.2f} seconds.")
    return computer


def solve_diffusion(
    u_data: pr.DataInterpolator,
    u0: df.Function,
    V: df.FunctionSpace,
    form: pr.TimedependentForm,
    coefficients: pr.CoefficientsDict,
    boundaries: list[pr.BoundaryData],
    time: pr.TimeKeeper,
    solver: pr.StationaryProblemSolver,
    storage: pr.FenicsStorage,
    computer: Optional[pr.BaseComputer] = None,
):
    computer = pr.set_optional(computer, pr.NullComputer)

    storage.write_function(u0, "total_concentration")
    computer.compute(time, u0)

    dirichlet_bcs = pr.process_dirichlet(V, boundaries)

    F = form(V, coefficients, boundaries, u0, time.dt)
    a = df.lhs(F)
    l = df.rhs(F)
    A = df.assemble(a)

    u = df.Function(V, name="total_concentration")
    tic = pytime.time()
    for ti in time:
        pr.print_progress(float(ti), time.endtime, rank=df.MPI.comm_world.rank)
        u_data.update(float(ti))
        b = df.assemble(l)
        solver.solve(u, A, b, dirichlet_bcs)
        computer.compute(ti, u)
        storage.write_checkpoint(u, "total_concentration", float(ti))
        u0.assign(u)

    storage.close()

    log("info", "Time loop finished.")
    toc = pytime.time()
    df.MPI.comm_world.barrier()
    log("info", f"Elapsed time in loop: {toc - tic:.2f} seconds.")
    return computer
