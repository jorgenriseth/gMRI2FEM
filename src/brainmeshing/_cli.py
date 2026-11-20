import click

from _cli import LazyGroup


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "meshgen": "brainmeshing.mesh_generation.meshgen",
        "process-surfaces": "brainmeshing.mesh_generation.process_surfaces",
        "extract-ventricles": "brainmeshing.mesh_generation.extract_ventricles",
        "refine_mesh": "brainmeshing.boundary_cell_refinement.refine_mesh_cli",
    },
)
def brainmeshing():
    pass
