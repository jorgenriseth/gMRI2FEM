import click
from gmri2fem._cli import mri, seg
from brainmeshing._cli import brainmeshing
from dti._cli import dti
from i2m._cli import i2m


@click.group()
def gmri2fem():
    pass


gmri2fem.add_command(mri)
gmri2fem.add_command(seg)
gmri2fem.add_command(dti)
gmri2fem.add_command(brainmeshing)
gmri2fem.add_command(i2m)
