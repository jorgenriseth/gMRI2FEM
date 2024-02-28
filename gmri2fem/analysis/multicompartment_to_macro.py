import logging
from pathlib import Path

import dolfin as df
from pantarei.fenicsstorage import FenicsStorage

from gmri2fem.models.multidiffusion_model import get_default_coefficients, print_progress

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
df.set_log_level(df.LogLevel.WARNING)


def multicompartment_fluid_concentrations_to_macro(
    input_data: Path,
    output_data: Path,
    compartments: list[str],
    porosities: dict[str, float] = None,
    store_name: str = None,
):
    if store_name is None:
        store_name = "multidiffusion_total"
    if porosities is None:
        porosities = get_default_coefficients()["porosity"]
    phi = [porosities[j] for j in compartments]

    logger.info(f"Reading 'multidiffusion' from {input_data}")
    inputfile = FenicsStorage(input_data, "r")
    timevec = inputfile.read_timevector("multidiffusion")
    u = inputfile.read_function("multidiffusion", idx=0)

    logger.info(f"Computing macroscopic concentration at t={timevec[0]}")
    c = df.Function(u.function_space().sub(0).collapse())
    c.vector()[:] = sum(
        [phi[j] * uj.vector() for j, uj in enumerate(u.split(deepcopy=True))]
    )
    outputfile = FenicsStorage(output_data, "w")
    outputfile.write_function(c, "multidiffusion", overwrite=True)

    for idx, ti in enumerate(timevec[1:]):
        print_progress(float(ti), timevec[-1], rank=df.MPI.comm_world.rank)
        u = inputfile.read_checkpoint(u, "multidiffusion", idx=idx + 1)
        c.vector()[:] = sum(
            [phi[j] * uj.vector() for j, uj in enumerate(u.split(deepcopy=True))]
        )
        outputfile.write_checkpoint(c, "multidiffusion_total", ti)
    inputfile.close()
    outputfile.close()

    file = FenicsStorage(outputfile.filepath, "r")
    file.to_xdmf("multidiffusion_total", "total")
    file.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)

    args = parser.parse_args()
    multicompartment_fluid_concentrations_to_macro(
        Path(args.input),
        Path(args.output),
        compartments=["ecs", "pvs"],
        porosities=None,
    )
