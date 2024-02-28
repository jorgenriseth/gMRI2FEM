#!/usr/bin/python3
"""
Usage:
python mesh_information.py 
    --paths [path/to/res1/mesh.hdf] [path/to/res2/mesh.hdf]
    --names res1 res2 # Optional, otherwise uses path as column index.
"""
import argparse
from pathlib import Path

import pantarei as pr

parser = argparse.ArgumentParser()
parser.add_argument("paths", nargs="+", type=Path)
parser.add_argument("--names", nargs="+", type=str)
args = parser.parse_args()


if hasattr(args, "names") and args.names is not None:
    if len(args.names) != len(args.paths):
        raise ValueError(
            f"paths-args and names-list must be same length, got {args.paths}, {args.names}"
        )
    indices = args.names
else:
    indices = list(map(str, args.paths))


meshes_info = {}
for index, hdf_file_path in zip(indices, args.paths):
    if not hdf_file_path.exists():
        raise FileNotFoundError(
            f"The specified HDF file does not exist: {hdf_file_path}"
        )
    mesh = pr.hdf2fenics(hdf_file_path, pack=True)
    num_cells = mesh.num_cells()
    num_vertices = mesh.num_vertices()
    min_diameter = 2 * mesh.rmin()
    max_diameter = 2 * mesh.rmax()
    mesh_quality = pr.utils.print_mesh_quality_info(mesh)
    meshes_info[hdf_file_path] = mesh_quality


indentsize = 4
indent = " " * indentsize
str_out = pr.utils.latex_header(index_header="Resolution")
for index, info in zip(indices, meshes_info.values()):
    str_out += 3 * indent
    str_out += pr.utils.latex_row(**info, decimals=3, index=index)
    str_out += "\n"
str_out += pr.utils.latex_footer()
print(str_out)
