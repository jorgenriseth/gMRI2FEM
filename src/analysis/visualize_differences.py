import argparse
from pathlib import Path

import dolfin as df
import numpy as np
import pantarei as pr

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timestamps", type=Path)
    parser.add_argument("--model_funcname", type=str, default="total_concentration")
    parser.add_argument("--reference_funcname", type=str, default="total_concentration")
    args = parser.parse_args()

    reference, model = args.reference, args.model

    # Load reference data.
    storage_reference = pr.FenicsStorage(reference, "r")
    storage_model = pr.FenicsStorage(model, "r")

    # Load data domain, to load into same value.
    u_ref = storage_reference.read_function(args.reference_funcname)
    V = u_ref.function_space()
    domain = V.mesh()
    if args.timestamps is not None:
        timevec = np.loadtxt(args.timestamps)
    else:
        timevec = storage_reference.read_timevector(args.reference_funcname)

    # Create new XDMF-file for storing difference
    with df.XDMFFile(df.MPI.comm_world, str(args.output)) as xdmf:
        for idx, ti in enumerate(timevec):
            pr.print_progress(float(ti), timevec[-1], rank=df.MPI.comm_world.rank)

            # Interpolate data from reference
            tvec = storage_reference.read_timevector(args.reference_funcname)
            bin = np.digitize(ti, tvec) - 1
            C = [
                storage_reference.read_function(
                    args.reference_funcname, domain=domain, idx=i
                )
                for i in range(tvec.size)[bin : bin + 2]
            ]
            interpolator_reference = pr.vectordata_interpolator(C, tvec[bin : bin + 2])
            # u_ref = interpolator_reference(ti)storage_reference.read_function(
            #     "total_concentration", domain=domain, idx=idx

            # Load storage data and interpolate to correct timepoint.
            tvec = storage_model.read_timevector(args.model_funcname)
            bin = np.digitize(ti, tvec) - 1
            C = [
                storage_model.read_function(args.model_funcname, domain=domain, idx=i)
                for i in range(tvec.size)[bin : bin + 2]
            ]
            interpolator_model = pr.vectordata_interpolator(C, tvec[bin : bin + 2])

            # Compute differences
            udiff = df.Function(V, name=f"{args.output.stem}")
            udiff.vector()[:] = interpolator_model(ti) - interpolator_reference(ti)
            print(type(V))
            exit()

            # And store it.
            xdmf.write(udiff, ti)

    storage_reference.close()
    storage_model.close()
