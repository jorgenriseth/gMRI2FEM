import click
from i2m.collect_mesh_data import collect
from i2m.concentrations_to_mesh import concentrations2mesh
from i2m.dti_data_to_mesh import dti2mesh
from i2m.mesh_segments import subdomains
from i2m.vtk_converter import hdf2vtk


@click.group()
def i2m():
    pass


i2m.add_command(collect)
i2m.add_command(concentrations2mesh)
i2m.add_command(dti2mesh)
i2m.add_command(subdomains)
i2m.add_command(hdf2vtk)
