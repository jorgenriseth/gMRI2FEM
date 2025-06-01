from pathlib import Path

import click
import numpy as np
import simple_mri as sm


@click.command()
@click.option("--input", type=Path, required=True)
@click.option("--reference", type=Path, required=True)
@click.option("--mask", type=Path, required=True)
@click.option("--output", type=Path, required=True)
def T1w_sigdiff(input: Path, reference: Path, mask: Path, output: Path):
    vol_mri = sm.load_mri(input, dtype=np.single)
    ref_mri = sm.load_mri(reference, dtype=np.single)
    vol = vol_mri.data
    ref = ref_mri.data

    mask_mri = sm.load_mri(mask, dtype=bool)
    signal_diff = compute_relative_signal_difference(vol, ref, mask_mri.data)
    signal_diff_mri = sm.SimpleMRI(signal_diff, affine=vol_mri.affine)
    sm.save_mri(signal_diff_mri, output, np.single)


def normalize_image(input: Path, refroi: Path, output: Path) -> Path:
    image = sm.load_mri(input, dtype=np.single)
    refroi_mri = sm.load_mri(refroi, dtype=bool)
    sm.assert_same_space(image, refroi_mri)
    normalized_data = image.data / np.median(image.data[refroi_mri.data])
    normalized_mri = sm.SimpleMRI(normalized_data, image.affine)
    sm.save_mri(normalized_mri, output, dtype=np.single)
    return output


def compute_relative_signal_difference(
    volume: np.ndarray, ref: np.ndarray, mask: np.ndarray
):
    nonzero_mask = mask * (ref > 0)
    signal_diff = np.nan * np.zeros_like(volume)
    signal_diff[nonzero_mask] = volume[nonzero_mask] / ref[nonzero_mask] - 1.0
    return signal_diff


@click.command()
@click.option("--input", type=Path, required=True)
@click.option("--refroi", type=Path, required=True)
@click.option("--output", type=Path, required=True)
def T1w_normalize(**kwargs):
    normalize_image(**kwargs)
