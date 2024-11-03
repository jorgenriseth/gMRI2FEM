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
import numpy as np
import nibabel

from nibabel.processing import resample_from_to


def find_valid_adjacent_tensor(data, i, j, k, max_iter):
    # Start at 1, since 0 is an invalid tensor
    for m in range(1, max_iter + 1):
        # Extract the adjacent data to voxel i, j, k
        # and compute the mean diffusivity.
        A = data[i - m : i + m + 1, j - m : j + m + 1, k - m : k + m + 1, :]
        A = A.reshape(-1, 9)
        MD = (A[:, 0] + A[:, 4] + A[:, 8]) / 3.0

        # If valid tensor is found:
        if MD.sum() > 0.0:
            # Find index of the median valid tensor, and return
            # corresponding tensor.
            index = (np.abs(MD - np.median(MD[MD > 0]))).argmin()
            return A[index]

    print("Failed to find valid tensor")
    return data[i, j, k]


def clean_dti_data(dti_file, mask_file, out_file, order=3, max_search=9):
    valid, mask, D = check_dti_data(dti_file, mask_file, order=order)
    # Zero out "invalid" tensor entries outside mask,
    # and extrapolate from valid neighbors
    D[~mask] = np.zeros(9)
    D[(~valid) * mask] = np.zeros(9)
    ii, jj, kk = np.where((~valid) * mask)
    for i, j, k in zip(ii, jj, kk):
        D[i, j, k, :] = find_valid_adjacent_tensor(D, i, j, k, max_search)

    # Create and save clean DTI image in T1 voxel space:
    mask_image = nibabel.load(mask_file)
    M1, M2, M3 = mask.shape
    shape = np.zeros((M1, M2, M3, 9))

    vox2ras = mask_image.header.get_vox2ras()
    dti_image = nibabel.nifti1.Nifti1Image(D, vox2ras)
    nibabel.nifti1.save(dti_image, out_file)

def compute_FA(eigvals: np.ndarray[float, (3,3)]):
    MD = (eigvals[:, 0] + eigvals[:, 1] + eigvals[:, 2]) / 3.0
    FA2 = (
        (3.0 / 2.0)
        * (
            (eigvals[:, 0] - MD) ** 2
            + (eigvals[:, 1] - MD) ** 2
            + (eigvals[:, 2] - MD) ** 2
        )
        / (eigvals[:, 0] ** 2 + eigvals[:, 1] ** 2 + eigvals[:, 2] ** 2)
    )
    FA = np.sqrt(FA2)
    return FA

def check_dti_data(dti_file, mask_file, order=0):
    # Load the DTI image data and mask:
    dti_image = nibabel.load(dti_file) dti_data = dti_image.get_fdata()

    mask_image = nibabel.load(mask_file)
    mask = mask_image.get_fdata().astype(bool)

    # Examine the differences in shape
    print("dti shape  ", dti_data.shape)
    print("mask shape ", mask.shape)
    M1, M2, M3 = mask.shape

    # Create an empty image as a helper for mapping
    # from DTI voxel space to T1 voxel space:
    shape = np.zeros((M1, M2, M3, 9))
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
    lmbdas, v = np.linalg.eigh(D)



    # Compute fractional anisotropy (FA)
    FA = compute_FA(lmbdas)

    # Define valid entries as those where all eigenvalues are
    # positive and FA is between 0 and 1
    positives = (lmbdas[:, 0] > 0) * (lmbdas[:, 1] > 0) * (lmbdas[:, 2] > 0)
    valid = positives * (FA < 1.0) * (FA > 0.0)
    valid = valid.reshape((M1, M2, M3))

    # Find all voxels with invalid tensors within the mask
    ii, jj, kk = np.where((~valid) * mask)
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dti", type=str)
    parser.add_argument("--mask", type=str)
    parser.add_argument("--out", type=str)
    parser.add_argument("--max_search", default=9, type=int)
    parser.add_argument("--order", default=3, type=int)

    Z = parser.parse_args()

    clean_dti_data(Z.dti, Z.mask, Z.out, Z.order, Z.max_search)
