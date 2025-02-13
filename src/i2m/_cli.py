import click
from _cli import LazyGroup


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "collect": "i2m.collect_mesh_data.collect",
        "concentrations2mesh": "i2m.concentrations_to_mesh.concentrations2mesh",
        "dti2mesh": "i2m.dti_data_to_mesh.dti2mesh",
        "subdomains": "i2m.mesh_segments.subdomains",
        "hdf2vtk": "i2m.vtk_converter.hdf2vtk",
        "vtk2mri": "i2m.vtk2mri.vtk2mri",
    },
)
def i2m():
    from i2m.collect_mesh_data import collect
    from i2m.concentrations_to_mesh import concentrations2mesh
    from i2m.dti_data_to_mesh import dti2mesh
    from i2m.mesh_segments import subdomains
    from i2m.vtk2mri import vtk2mri
    from i2m.vtk_converter import hdf2vtk

    i2m.add_command(collect)
    i2m.add_command(concentrations2mesh)
    i2m.add_command(dti2mesh)
    i2m.add_command(subdomains)
    i2m.add_command(hdf2vtk)
    i2m.add_command(vtk2mri)
