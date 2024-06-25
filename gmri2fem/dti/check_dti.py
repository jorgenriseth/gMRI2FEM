"""
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
"""

import argparse
import numpy
import nibabel

from nibabel.processing import resample_from_to

numpy.seterr(divide="ignore", invalid="ignore")


def check_dti_data(dti_file, mask_file, order=0):
    # Load the DTI image data and mask:
    dti_image = nibabel.load(dti_file)
    dti_data = dti_image.get_fdata()

    mask_image = nibabel.load(mask_file)
    mask = mask_image.get_fdata().astype(bool)

    # Examine the differences in shape
    print("dti shape  ", dti_data.shape)
    print("mask shape ", mask.shape)
    M1, M2, M3 = mask.shape

    # Create an empty image as a helper for mapping
    # from DTI voxel space to T1 voxel space:
    shape = numpy.zeros((M1, M2, M3, 9))
    vox2ras = mask_image.header.get_vox2ras()
    Nii = nibabel.nifti1.Nifti1Image
    helper = Nii(shape, vox2ras)

    # Resample the DTI data in the T1 voxel space:
    image = resample_from_to(dti_image, helper, order=order)
    D = image.get_fdata()
    print("resampled image shape ", D.shape)

    # Reshape D from M1 x M2 x M3 x 9 into a N x 3 x 3:
    D = D.reshape(-1, 3, 3)

    # Compute eigenvalues and eigenvectors
    lmbdas, v = numpy.linalg.eigh(D)

    def compute_FA(lmbdas):
        MD = (lmbdas[:, 0] + lmbdas[:, 1] + lmbdas[:, 2]) / 3.0
        FA2 = (
            (3.0 / 2.0)
            * (
                (lmbdas[:, 0] - MD) ** 2
                + (lmbdas[:, 1] - MD) ** 2
                + (lmbdas[:, 2] - MD) ** 2
            )
            / (lmbdas[:, 0] ** 2 + lmbdas[:, 1] ** 2 + lmbdas[:, 2] ** 2)
        )
        FA = numpy.sqrt(FA2)
        return FA

    # Compute fractional anisotropy (FA)
    FA = compute_FA(lmbdas)

    # Define valid entries as those where all eigenvalues are
    # positive and FA is between 0 and 1
    positives = (lmbdas[:, 0] > 0) * (lmbdas[:, 1] > 0) * (lmbdas[:, 2] > 0)
    valid = positives * (FA < 1.0) * (FA > 0.0)
    valid = valid.reshape((M1, M2, M3))

    # Find all voxels with invalid tensors within the mask
    ii, jj, kk = numpy.where((~valid) * mask)
    print("Number of invalid tensor voxels within the mask ROI: ", len(ii))

    # Reshape D from N x 3 x 3 to M1 x M2 x M3 x 9
    D = D.reshape((M1, M2, M3, 9))

    return valid, mask, D


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dti", type=str)
    parser.add_argument("--mask", type=str)
    parser.add_argument("--order", default=0, type=int)

    Z = parser.parse_args()

    check_dti_data(Z.dti, Z.mask, Z.order)
