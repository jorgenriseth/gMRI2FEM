import click
from _cli import LazyGroup


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "mesh-generation": "brainmeshing.mesh_generation.meshgen",
        "process-surfaces": "brainmeshing.mesh_generation.process_surfaces",
    },
)
def dti():
    from dti.reslice_dti import reslice_dti
    from dti.clean_dti_data import clean
    from dti.utils import create_eddy_index_file, create_topup_params_file
