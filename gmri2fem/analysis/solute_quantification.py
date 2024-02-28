import logging
from pathlib import Path
from typing import Callable

import dolfin as df
from pantarei.computers import BaseComputer
from pantarei.fenicsstorage import FenicsStorage
from tqdm import tqdm

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
df.set_log_level(df.LogLevel.WARNING)


def solute_quantification(
    inputpath: Path,
    funcname: str,
    computer_factory: Callable[[df.Measure], BaseComputer],
) -> BaseComputer:
    inputfile = FenicsStorage(inputpath, "r")
    timevec = inputfile.read_timevector(funcname)
    u = inputfile.read_function(funcname, idx=0)
    domain = u.function_space().mesh()
    ds = df.Measure("ds", domain=domain, subdomain_data=domain.boundaries)
    dx = df.Measure("dx", domain=domain, subdomain_data=domain.subdomains)
    computer = computer_factory(dx, ds)
    computer.init_from_vector(timevec)

    for idx, ti in tqdm(enumerate(timevec)):
        u = inputfile.read_checkpoint(u, funcname, idx)
        computer.compute_from_index(idx, u)
    inputfile.close()
    return computer


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--funcname", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()

    create_computer = lambda dx, ds: BaseComputer(
        {
            "whole-brain": lambda u: df.assemble(u * dx),
            "gray-matter": lambda u: df.assemble(u * dx(1)),
            "white-matter": lambda u: df.assemble(u * dx(2)),
            "pial-surf": lambda u: df.assemble(u * ds(4)),
            "inferior-surf": lambda u: df.assemble(u * ds(5)),
            "ventricle-surf": lambda u: df.assemble(u * ds(8)),
        }
    )
    computer = solute_quantification(Path(args.input), args.funcname, create_computer)

    import pandas as pd

    dframe = pd.DataFrame.from_dict(computer.values)
    dframe.to_csv(Path(args.output))
