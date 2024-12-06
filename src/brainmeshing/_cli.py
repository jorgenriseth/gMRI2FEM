import click

from brainmeshing.gray_matter_surfaces import main as gray_matter_surfaces
from brainmeshing.mesh_generation import meshgen as mesh_generation

from brainmeshing.ventricles import main as extract_ventricles
from brainmeshing.white_gray_separation import main as white_gray_separation
from brainmeshing.white_matter_surfaces import main as white_matter_surfaces


@click.group()
def brainmeshing():
    pass


brainmeshing.add_command(extract_ventricles)
brainmeshing.add_command(white_gray_separation)
brainmeshing.add_command(white_matter_surfaces)
brainmeshing.add_command(gray_matter_surfaces)
brainmeshing.add_command(mesh_generation)
