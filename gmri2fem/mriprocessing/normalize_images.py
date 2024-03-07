from pathlib import Path

import nibabel
import numpy as np


def normalize_image(image_path: Path, refroi_path: Path, outputpath: Path) -> Path:
    image = nibabel.freesurfer.mghformat.load(image_path)
    refroi = nibabel.freesurfer.mghformat.load(refroi_path)

    assert np.allclose(
        refroi.affine, image.affine
    ), "Poor match between reference and image-transform."
    image_data = image.get_fdata(dtype=np.float32)
    ref_mask = refroi.get_fdata().astype(bool)

    normalized_image_data = image_data / np.median(image_data[ref_mask])
    normalized_image = nibabel.freesurfer.mghformat.MGHImage(
        normalized_image_data, image.affine
    )
    nibabel.freesurfer.mghformat.save(normalized_image, outputpath)
    return outputpath


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--image", type=Path, required=True)
    parser.add_argument("--refroi", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    normalize_image(args.image, args.refroi, args.output)
