import click

from _cli import LazyGroup


# @click.group(
#     cls=LazyGroup,
#     lazy_subcommands={
#         "mesh-generation": "brainmeshing.mesh_generation.meshgen",
#         "process-surfaces": "brainmeshing.mesh_generation.process_surfaces",
#         "extract-ventricles": "brainmeshing.mesh_generation.extract_ventricles",
#     },
# )
@click.group()
def brainmeshing():
    pass


from brainmeshing.mesh_generation import extract_ventricles, meshgen, process_surfaces
from brainmeshing.boundary_cell_refinement import refine_mesh_cli

brainmeshing.add_command(meshgen)
brainmeshing.add_command(process_surfaces)
brainmeshing.add_command(extract_ventricles)
brainmeshing.add_command(refine_mesh_cli)
