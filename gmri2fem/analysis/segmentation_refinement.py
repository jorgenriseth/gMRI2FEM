import itertools
from pathlib import Path

import nibabel.nifti1 as nifti1
import nibabel.freesurfer.mghformat as mghformat
import numpy as np
import scipy
import skimage

from gmri2fem.utils import apply_affine


def seg_upsampling(
    reference: Path,
    segmentation: Path,
):
    seg_mgz = mghformat.load(segmentation)
    seg = seg_mgz.get_fdata(dtype=np.single).astype(np.int32)

    T1map_nii = nifti1.load(reference)
    T1map = T1map_nii.get_fdata(dtype=np.single)

    shape_in = seg.shape
    shape_out = T1map.shape

    upsampled_inds = np.fromiter(
        itertools.product(*(np.arange(ni) for ni in shape_out)),
        dtype=np.dtype((int, 3)),
    )

    seg_affine = seg_mgz.affine
    T1map_affine = T1map_nii.affine
    seg_inds = apply_affine(
        np.linalg.inv(seg_affine), apply_affine(T1map_affine, upsampled_inds)
    )
    seg_inds = np.rint(seg_inds).astype(np.int32)

    # The two images does not necessarily share field of view.
    # Remove voxels which are not located within the segmentation fov.
    valid_index_mask = (seg_inds > 0).all(axis=1) * (seg_inds < shape_in).all(axis=1)
    upsampled_inds = upsampled_inds[valid_index_mask]
    seg_inds = seg_inds[valid_index_mask]

    I_in, J_in, K_in = seg_inds.T
    I_out, J_out, K_out = upsampled_inds.T

    seg_upsampled = np.zeros(T1map.shape)
    seg_upsampled[I_out, J_out, K_out] = seg[I_in, J_in, K_in]
    return nifti1.Nifti1Image(seg_upsampled, T1map_affine)


def csf_segmentation(
    seg_upsampled_mri: nifti1.Nifti1Image,
    csf_mask_mri: nifti1.Nifti1Image,
) -> nifti1.Nifti1Image:
    seg_upsampled = seg_upsampled_mri.get_fdata(dtype=np.single).astype(np.int32)
    I, J, K = np.where(seg_upsampled != 0)
    inds = np.array([I, J, K]).T
    interp = scipy.interpolate.NearestNDInterpolator(inds, seg_upsampled[I, J, K])

    csf_mask = csf_mask_mri.get_fdata(dtype=np.single).astype(bool)
    i, j, k = np.where(csf_mask)

    csf_seg = np.zeros_like(seg_upsampled)
    csf_seg[i, j, k] = interp(i, j, k)
    return nifti1.Nifti1Image(csf_seg, csf_mask_mri.affine, header=csf_mask_mri.header)


def segmentation_refinement(
    upsampled_segmentation: nifti1.Nifti1Image,
    csf_segmentation: nifti1.Nifti1Image,
    closing_radius: int = 5,
) -> nifti1.Nifti1Image:
    seg_upsampled = upsampled_segmentation.get_fdata(dtype=np.single).astype(np.int32)

    combined_segmentation = seg_upsampled.copy()
    combined_segmentation = skimage.segmentation.expand_labels(
        combined_segmentation, distance=3
    )
    csf_seg = csf_segmentation.get_fdata(dtype=np.single).astype(np.int32)
    csf_mask = (csf_seg != 0).astype(bool)
    combined_segmentation[csf_mask] = -csf_seg[csf_mask]

    radius = closing_radius
    combined_mask = csf_mask + (seg_upsampled != 0)
    combined_mask = skimage.morphology.closing(
        combined_mask,
        footprint=np.ones([1 + radius * 2] * combined_mask.ndim),
    )
    combined_segmentation[~combined_mask] = 0
    aseg_new = np.where(combined_segmentation > 0, combined_segmentation, 0)
    return nifti1.Nifti1Image(aseg_new, upsampled_segmentation.affine)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--fs_seg", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--csfmask", type=Path, required=True)
    parser.add_argument("--output_seg", type=Path, required=True)
    parser.add_argument("--output_csfseg", type=Path, required=True)
    args = parser.parse_args()

    upsampled_seg = seg_upsampling(args.reference, args.fs_seg)
    csf_mask = nifti1.load(args.csfmask)
    csf_seg = csf_segmentation(upsampled_seg, csf_mask)
    nifti1.save(csf_seg, args.output_seg)

    refined_seg = segmentation_refinement(upsampled_seg, csf_seg)
    nifti1.save(refined_seg, args.output_csfseg)
