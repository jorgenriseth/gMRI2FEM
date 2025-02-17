import click
from _cli import LazyGroup


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "mesh-generation": "brainmeshing.mesh_generation.meshgen",
        "process-surfaces": "brainmeshing.mesh_generation.process_surfaces",
    },
)
def brainmeshing():
    pass
