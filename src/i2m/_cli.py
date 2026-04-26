import click

from _cli import LazyGroup


@click.group(
    cls=LazyGroup,
    lazy_subcommands={
        "collect": "i2m.collect_mesh_data.collect",
        "concentrations2mesh": "i2m.concentrations_to_mesh.concentrations2mesh",
        "tissue-concentrations": "i2m.voxel_center_minimization.map_evaluation_data_to_mesh",
        "boundary-concentrations": "i2m.concentrations_to_mesh.map_boundary_concentrations_cli",
        "dti2mesh": "i2m.dti_data_to_mesh.dti2mesh",
        "subdomains": "i2m.mesh_segments.subdomains",
        "hdf2vtk": "i2m.vtk_converter.hdf2vtk",
        "vtk2mri": "i2m.vtk2mri.vtk2mri",
        "evaluation_data": "i2m.evaluation_data.create_mesh_evaluation_data",
    },
)
def i2m():
    pass
