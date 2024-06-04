from pathlib import Path

import matplotlib.pyplot as plt
import nibabel
import numpy as np
import scipy
import skimage


def grow_refroi(
    image_data: np.ndarray, mask: np.ndarray, median_rel_tolerance: float = 0.10
):
    offsets = np.abs(
        image_data[:, mask] - np.median(image_data[:, mask], axis=1).reshape(-1, 1)
    )
    ind = np.array(np.where(mask)).T
    seed_voxels = ind[np.argmin(offsets, axis=1)]
    tolerances = median_rel_tolerance * np.median(image_data[:, mask], axis=1)
    refrois = np.array(
        [
            skimage.segmentation.flood(
                image_data[i],
                tuple(seed_voxels[i]),
                tolerance=tolerances[i],
                # connectivity=1,
            )
            for i in range(len(image_data))
        ]
    )
    refroi = np.all(refrois, axis=0)
    return refroi


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--references", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rel_tolerance", type=float, default=0.1)
    args = parser.parse_args()

    seedimage = nibabel.freesurfer.mghformat.load(args.seed)
    data = np.array(
        [
            nibabel.freesurfer.mghformat.load(im).get_fdata(dtype=np.float32)
            for im in args.references
        ]
    )
    maskseed = np.array(seedimage.get_fdata()) == 1.0
    refroi = grow_refroi(data, maskseed, median_rel_tolerance=args.rel_tolerance)
    im_out = nibabel.freesurfer.mghformat.MGHImage(
        refroi, seedimage.affine, seedimage.header
    )
    nibabel.freesurfer.mghformat.save(im_out, args.output)
