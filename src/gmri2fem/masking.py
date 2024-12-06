from pathlib import Path

import click
import numpy as np
import skimage
from simple_mri import load_mri, save_mri, SimpleMRI, assert_same_space

from gmri2fem.utils import largest_island


def intracranial_mask(csf_mask: SimpleMRI, segmentation: SimpleMRI):
    assert_same_space(csf_mask, segmentation)
    combined_mask = csf_mask.data + segmentation.data.astype(bool)
    background_mask = largest_island(~combined_mask, connectivity=1)
    opened = skimage.morphology.binary_opening(
        background_mask, skimage.morphology.ball(3)
    )
    return SimpleMRI(data=~opened, affine=segmentation.affine)


@click.command()
@click.option("--csfmask", type=Path, required=True)
@click.option("--brain_seg", "brain_mask", type=Path, required=True)
@click.option("--output", type=Path, required=True)
def mask_intracranial(csfmask: Path, brain_mask: Path, output):
    csf = load_mri(csfmask, dtype=bool)
    brain = load_mri(brain_mask, dtype=bool)
    intracranial_mask_mri = intracranial_mask(csf, brain)
    save_mri(intracranial_mask_mri, output, dtype=np.uint8)


@click.command()
@click.option("--input", type=Path, required=True)
@click.option("--output", type=Path, required=True)
@click.option("--connectivity", type=int, default=2)
@click.option("--li", type=bool, is_flag=True, default=False)
def mask_csf(input: Path, output: Path, connectivity: int = 2, li: bool = False):
    mri = load_mri(input, np.single)
    mask = create_csf_mask(mri.data, connectivity, li)

    assert np.max(mask) > 0, "Masking failed, no voxels in mask"
    save_mri(SimpleMRI(mask, mri.affine), output, np.uint8)


def create_csf_mask(
    vol: np.ndarray,
    connectivity: int = 2,
    use_li: bool = False,
) -> np.ndarray:
    if use_li:
        thresh = skimage.filters.threshold_li(vol)
        binary = vol > thresh
        binary = largest_island(binary, connectivity=connectivity)
    else:
        (hist, bins) = np.histogram(
            vol[(vol > 0) * (vol < np.quantile(vol, 0.999))], bins=512
        )
        thresh = skimage.filters.threshold_yen(hist=(hist, bins))
        binary = vol > thresh
        binary = largest_island(binary, connectivity=connectivity)
    return binary
