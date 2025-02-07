import click
from _cli import LazyGroup


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "mesh_generation": "brainmeshing.mesh_generation.meshgen",
        "process_surfaces": "brainmeshing.mesh_generation.process_surfaces",
    },
)
def brainmeshing():
    from brainmeshing.mesh_generation import process_surfaces
    from brainmeshing.mesh_generation import meshgen as mesh_generation
