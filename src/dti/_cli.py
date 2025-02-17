import click
from _cli import LazyGroup


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "reslice-dti": "dti.reslice_dti.reslice_dti",
        "clean": "dti.clean_dti_data.clean",
        "eddy-index": "dti.utils.create_eddy_index_file",
        "topup-params": "dti.utils.create_topup_params_file",
    },
)
def dti():
    pass
