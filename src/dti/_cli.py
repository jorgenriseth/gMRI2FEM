from pathlib import Path
from typing import Required
import click
from _cli import LazyGroup


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "reslice-dti": "dti.reslice_dti.reslice_dti",
        "clean": "dti.clean_dti_data.clean",
        "eddy-index": "dti.eddy.create_eddy_index_file",
        "eddy": "dti._cli.eddy_cli",
        "topup": "dti._cli.topup_cli",
    },
)
def dti():
    pass


@click.command()
@click.option("--input", "-i", "dti", type=Path, required=True)
@click.option("--topup", "-t", "topup_b0", type=Path, required=True)
@click.option("--outputdir", "-o", type=Path, required=True)
@click.option("--tmppath", type=Path)
def topup_cli(**kwargs):
    from dti.topup import topup

    topup(**kwargs)


@click.command()
@click.option("--input", "-i", "dti", type=Path, required=True)
@click.option("--topup", "-t", "topup_b0_mean", type=Path, required=True)
@click.option("--acq_params", "-a", type=Path, required=True)
@click.option("--output", "-o", type=Path, required=True)
@click.option("--multiband_factor", "--mb", type=int, default=1)
@click.option("--tmppath", type=Path)
def eddy_cli(**kwargs):
    from dti.eddy import eddy_correct

    eddy_correct(**kwargs)
