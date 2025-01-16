import click

from brainmeshing.mesh_generation import meshgen as mesh_generation


@click.group()
def brainmeshing():
    pass


brainmeshing.add_command(mesh_generation)
