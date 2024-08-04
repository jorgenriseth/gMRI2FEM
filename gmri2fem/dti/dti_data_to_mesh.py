"""
Adjusted version of script from original MRI2FEM-book according to below license.


'''
Copyright (c) 2020 Kent-Andre Mardal, Marie E. Rognes, Travis B. Thompson, Lars Magnus Valnes

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''
"""
import argparse
import numpy
import nibabel  # importing nibabel after dolfin will cause error
from nibabel.affines import apply_affine
import dolfin as df
from dolfin import *

import pantarei as pr


def adjusting_mean_diffusivity(Dvector, subdomains, tags_with_limits):
    # Computes the mean diffusivity for each degree of freedom
    MD = (Dvector[:, 0] + Dvector[:, 4] + Dvector[:, 8]) / 3.0

    # Reads the tag and the minimum and maximum mean diffusivity limit
    # for that tag subdomain.
    for tag, mn, mx in tags_with_limits:
        # If the minimum or maximum mean diffusivity limit is set to zero,
        # then the limit is considered void.
        usr_max = float(mx) if mx != 0 else numpy.inf
        usr_min = float(mn) if mn != 0 else -numpy.inf

        # creates a mask for all degrees of freesom that are within the
        # subdomain with the tag and is above the maximum limit or
        # below the minimum limit.
        max_mask = (subdomains.array() == tag) * (MD > usr_max)
        min_mask = (subdomains.array() == tag) * (MD < usr_min)

        # Sets values that are either above or below limits to the closest limit.
        Dvector[max_mask] = usr_max * numpy.divide(
            Dvector[max_mask], MD[max_mask, numpy.newaxis]
        )
        Dvector[min_mask] = usr_min * numpy.divide(
            Dvector[min_mask], MD[min_mask, numpy.newaxis]
        )


def dti_data_to_mesh(meshfile, dti, outfile, label=None):
    # Read the mesh from file. The mesh coordinates define
    # the Surface RAS space.
    mesh = Mesh()
    hdf = HDF5File(mesh.mpi_comm(), meshfile, "r")
    hdf.read(mesh, "domain/mesh", False)

    # Read subdomains and boundary markers to write output.
    d = mesh.topology().dim()
    subdomains = MeshFunction("size_t", mesh, d)
    hdf.read(subdomains, "domain/subdomains")
    boundaries = MeshFunction("size_t", mesh, d - 1)
    hdf.read(boundaries, "domain/boundaries")
    hdf.close()

    # Read in DTI data in T1 voxel space
    # (see previous function, e.g. ernie-dti-clean.mgz)
    dti_image = nibabel.load(dti)
    dti_data = dti_image.get_fdata()

    # Transformation to voxel space from mesh coordinates
    # vox2ras = dti_image.header.get_vox2ras()
    vox2ras = dti_image.affine
    ras2vox = numpy.linalg.inv(vox2ras)

    # Create a FEniCS tensor field:
    DG09 = TensorFunctionSpace(mesh, "DG", 0)
    D = Function(DG09)

    # Get the coordinates xyz of each degree of freedom
    DG0 = FunctionSpace(mesh, "DG", 0)
    imap = DG0.dofmap().index_map()
    num_dofs_local = imap.local_range()[1] - imap.local_range()[0]
    xyz = DG0.tabulate_dof_coordinates()
    xyz = xyz.reshape((num_dofs_local, -1))

    # Convert to voxel space and round off to find
    # voxel indices
    ijk = apply_affine(ras2vox, xyz).T
    i, j, k = numpy.rint(ijk).astype("int")

    # Create a matrix from the DTI representation
    D1 = dti_data[i, j, k]
    print(D1.shape)

    # Further manipulate data (described better later)
    if label:
        adjusting_mean_diffusivity(D1, subdomains, label)

    # Assign the output to the tensor function
    D.vector()[:] = D1.reshape(-1)

    # Compute other functions
    md = 1.0 / 3.0 * tr(D)
    MD = project(md, DG0, solver_type="cg", preconditioner_type="amg")
    fa = sqrt((3.0 / 2) * inner(dev(D), dev(D)) / inner(D, D))
    FA = project(fa, DG0)

    # Now store everything to a new file - ready for use!
    domain = pr.Domain(mesh, subdomains, boundaries)
    hdf = HDF5File(mesh.mpi_comm(), outfile, "w")
    pr.write_domain(hdf, domain)
    pr.write_function(hdf, D, "DTI")
    pr.write_function(hdf, MD, "MD")
    pr.write_function(hdf, FA, "FA")
    hdf.close()


if __name__ == "__main__":
    # Create parser for command line arguments for files and filenames
    # TODO: Write about arguments or simplify in the book!
    parser = argparse.ArgumentParser()
    parser.add_argument("--dti", type=str, default="")
    parser.add_argument("--mesh", type=str)
    parser.add_argument("--out", type=str)
    parser.add_argument(
        "--label",
        type=float,
        action="append",
        nargs=3,
        help="--label TAG MIN MAX. The value of zero is considered void.",
    )
    Z = parser.parse_args()

    dti_data_to_mesh(Z.mesh, Z.dti, Z.out, Z.label)

    # In the book, just add line with:
    # dti_data_to_mesh("ernie.h5", "ernie-dti-clean.mgz", "ernie-dti.h5")
