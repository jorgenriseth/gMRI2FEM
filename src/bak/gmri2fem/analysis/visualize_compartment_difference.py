import argparse
from pathlib import Path

import dolfin as df
import numpy as np
import pantarei as pr


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--funcname")
    parser.add_argument("--fname1", type=Path, required=True)
    parser.add_argument("--fname2", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    inputfile = args.input
    # Load reference data.
    storage = pr.FenicsStorage(inputfile, "r")

    # Load data domain, to load into same value.
    u_ref = storage.read_function(
        args.funcname
    )
    W = u_ref.function_space()
    domain = W.mesh()
    timevec = storage.read_timevector(
        args.funcname
    )
    u_diff = df.Function(W.sub(0).collapse(), name=f"{args.fname1}-{args.fname2}")


    # Create new XDMF-file for storing difference
    with df.XDMFFile(df.MPI.comm_world, str(args.output)) as xdmf:
        for idx, ti in enumerate(timevec):
            pr.print_progress(float(ti), timevec[-1], rank=df.MPI.comm_world.rank)

            # Interpolate data from reference
            u = storage.read_function(args.funcname, domain, idx)
            u1, u2 = u.split(deepcopy=True)
            
            # Compute differences
            u_diff.vector()[:] = u1.vector()[:] - u2.vector()[:]
            
            # And store it.
            xdmf.write(u_diff, ti)

    storage.close()
