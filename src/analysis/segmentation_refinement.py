import itertools
from pathlib import Path

import nibabel.freesurfer.mghformat as mghformat
import nibabel.nifti1 as nifti1
import numpy as np
import scipy
import skimage

from gmri2fem.utils import apply_affine

MASK_DTYPE = np.uint8
SEG_DTYPE = np.int16
DATA_DTYPE = np.single


def seg_upsampling(
    reference: Path,
    segmentation: Path,
):
    seg_mgz = mghformat.load(segmentation)
    seg = np.asanyarray(seg_mgz.dataobj, dtype=SEG_DTYPE)

    reference_nii = nifti1.load(reference)
    shape_in = seg.shape
    shape_out = reference_nii.shape

    upsampled_inds = np.fromiter(
        itertools.product(*(np.arange(ni) for ni in shape_out)),
        dtype=np.dtype((int, 3)),
    )

    seg_affine = seg_mgz.affine
    reference_affine = reference_nii.affine
    seg_inds = apply_affine(
        np.linalg.inv(seg_affine), apply_affine(reference_affine, upsampled_inds)
    )
    seg_inds = np.rint(seg_inds).astype(SEG_DTYPE)

    # The two images does not necessarily share field of view.
    # Remove voxels which are not located within the segmentation fov.
    valid_index_mask = (seg_inds > 0).all(axis=1) * (seg_inds < shape_in).all(axis=1)
    upsampled_inds = upsampled_inds[valid_index_mask]
    seg_inds = seg_inds[valid_index_mask]

    I_in, J_in, K_in = seg_inds.T
    I_out, J_out, K_out = upsampled_inds.T

    seg_upsampled = np.zeros(shape_out, dtype=SEG_DTYPE)
    seg_upsampled[I_out, J_out, K_out] = seg[I_in, J_in, K_in]
    return nifti1.Nifti1Image(seg_upsampled, reference_affine)


def segment_csf(
    seg_upsampled_mri: nifti1.Nifti1Image,
    csf_mask_mri: nifti1.Nifti1Image,
) -> nifti1.Nifti1Image:
    seg_upsampled = np.asanyarray(seg_upsampled_mri.dataobj, dtype=SEG_DTYPE)
    I, J, K = np.where(seg_upsampled != 0)
    inds = np.array([I, J, K]).T
    interp = scipy.interpolate.NearestNDInterpolator(inds, seg_upsampled[I, J, K])

    csf_mask = np.asanyarray(csf_mask_mri.dataobj, dtype=bool)
    i, j, k = np.where(csf_mask)

    csf_seg = np.zeros_like(seg_upsampled)
    csf_seg[i, j, k] = interp(i, j, k)
    return nifti1.Nifti1Image(csf_seg, csf_mask_mri.affine)


def segmentation_refinement(
    upsampled_segmentation: nifti1.Nifti1Image,
    csf_segmentation: nifti1.Nifti1Image,
    closing_radius: int = 5,
) -> nifti1.Nifti1Image:
    seg_upsampled = np.asanyarray(upsampled_segmentation.dataobj, dtype=SEG_DTYPE)

    combined_segmentation = seg_upsampled.copy()
    combined_segmentation = skimage.segmentation.expand_labels(
        combined_segmentation, distance=3
    )
    csf_seg = np.asanyarray(csf_segmentation.dataobj, dtype=SEG_DTYPE)
    csf_mask = csf_seg != 0
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
    parser.add_argument("--intracranial_mask", type=Path)
    args = parser.parse_args()

    csf_mask = nifti1.load(args.csfmask)
    upsampled_seg = seg_upsampling(args.reference, args.fs_seg)
    csf_seg = segment_csf(upsampled_seg, csf_mask)
    nifti1.save(csf_seg, args.output_csfseg)

    refined_seg = segmentation_refinement(upsampled_seg, csf_seg)
    nifti1.save(refined_seg, args.output_seg)

    if args.intracranial_mask is not None:
        csf_mask_data = csf_mask.get_fdata().astype(bool)
        refined_seg_data = refined_seg.get_fdata().astype(bool)
        intracranial = csf_mask_data + (refined_seg_data != 0)
        intracranial_nii = nifti1.Nifti1Image(
            intracranial.astype(np.single), refined_seg.affine
        )
        nifti1.save(intracranial_nii, args.intracranial_mask)
