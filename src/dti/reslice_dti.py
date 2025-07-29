import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import click
import numpy as np
from simple_mri import SimpleMRI, load_mri, save_mri

from dti.utils import mri_number_of_frames
from gmri2fem.reslice_4d import reslice_4d


def construct_tensor_from_eigs(
    dti_folder: Path, prefix_pattern: str, suffix: str
) -> SimpleMRI:
    eigvec = load_mri(dti_folder / f"{prefix_pattern}_V1{suffix}.nii.gz", np.single)
    spatial_shape = eigvec.data.shape[:3]
    B = np.zeros((*spatial_shape, 3, 3), dtype=eigvec.data.dtype)
    L = np.zeros_like(B)
    for i in range(1, 4):
        eigvec = load_mri(
            dti_folder / f"{prefix_pattern}_V{i}{suffix}.nii.gz", np.single
        )
        eigval = load_mri(
            dti_folder / f"{prefix_pattern}_L{i}{suffix}.nii.gz", np.single
        )
        B[..., :, i - 1] = eigvec.data
        L[..., i - 1, i - 1] = eigval.data
    valid_mask = np.linalg.det(L) != 0
    Binv = np.zeros_like(B)
    Binv[valid_mask] = np.linalg.inv(B[valid_mask])
    return SimpleMRI((B @ L @ Binv).reshape(*spatial_shape, 9), eigvec.affine)


def construct_tensor_from_vector_array(tensor: SimpleMRI) -> SimpleMRI:
    spatial_shape = tensor.data.shape[:3]
    Dt = np.zeros((*spatial_shape, 3, 3), dtype=tensor.data.dtype)
    if tensor.data.shape[-1] == 6:
        Dt[..., 0, :] = Dt[..., :, 0] = tensor.data[..., :3]
        Dt[..., 1, 1:] = Dt[..., 1:, 1] = tensor.data[..., 3:5]
        Dt[..., 2, 2] = tensor.data[..., 5]
    elif tensor.data.shape[-1] == 9:
        Dt[...] = tensor.data.reshape(*spatial_shape, 3, 3).copy()
    return SimpleMRI(Dt, tensor.affine)


@click.command()
@click.option("-input_folder", type=Path, required=True)
@click.option("-prefix_pattern", type=Path, required=True)
@click.option("-output", type=Path, required=True)
@click.option("--from_tensor", default=False)
def construct_and_save_tensor(
    input_folder: Path, prefix_pattern: str, output: Path, from_tensor: bool
):
    if from_tensor:
        tensor_in = load_mri(
            input_folder / f"{prefix_pattern}_tensor.nii.gz", dtype=np.single
        )
        tensor_out = construct_tensor_from_vector_array(tensor_in)
    else:
        tensor_out = construct_tensor_from_eigs(input_folder, prefix_pattern, "")
    save_mri(tensor_out, output, dtype=tensor_out.data.dtype)


@click.command()
@click.option("--fixed", type=Path, required=True)
@click.option("--dtidir", type=Path, required=True)
@click.option("--prefix_pattern", type=str, required=True)
@click.option("--outdir", type=Path)
@click.option("--transform", type=Path)
@click.option("--threads", type=int, default=1)
@click.option("--suffix", type=str, default="")
@click.option("--interp_mode", type=str, default="NN")
@click.option("--greedyargs", type=str, default="")
def reslice_dti(
    fixed: Path,
    dtidir: Path,
    prefix_pattern: str,
    outdir: Path,
    transform: Path,
    threads: int,
    out_pattern: Optional[str] = None,
    suffix: str = None,
    interp_mode: str = "NN",
    greedyargs: str = None,
):
    if out_pattern is None:
        out_pattern = prefix_pattern

    if greedyargs is None:
        greedyargs = ""

    for c in ["FA", "MD", *[f"L{i}" for i in range(1, 4)]]:
        inpath = dtidir / f"{prefix_pattern}_{c}.nii.gz"
        outpath = outdir / f"{out_pattern}_{c}{suffix}.nii.gz"
        reslice_4d(
            inpath,
            fixed,
            outpath,
            transform,
            threads,
            interp_mode=interp_mode,
            greedyargs=greedyargs,
        )

    with tempfile.TemporaryDirectory(prefix=out_pattern) as tmpdir:
        tmppath = Path(tmpdir)
        for Vi in [f"V{i}" for i in range(1, 4)]:
            inpath = dtidir / f"{prefix_pattern}_{Vi}.nii.gz"
            resliced = tmppath / f"{out_pattern}_{Vi}{suffix}.nii.gz"
            outpath = outdir / resliced.name
            reslice_4d(
                inpath,
                fixed,
                resliced,
                transform,
                threads,
                interp_mode=interp_mode,
                greedyargs=greedyargs,
            )
            resliced_mri = load_mri(resliced, dtype=np.single)
            norms = np.linalg.norm(resliced_mri.data, axis=-1, ord=2)
            resliced_mri.data[norms > 0] /= norms[norms > 0, np.newaxis]
            save_mri(resliced_mri, outpath, dtype=np.single)

    resliced_tensor = construct_tensor_from_eigs(outdir, out_pattern, suffix)
    save_mri(
        resliced_tensor,
        outdir / f"{out_pattern}_tensor{suffix}.nii.gz",
        dtype=np.single,
    )


@click.command("eddy-index")
@click.option("--input", type=Path, required=True)
@click.option("--output", type=Path, required=True)
def create_eddy_index_file(input: Path, output: Path):
    nframes = mri_number_of_frames(input)
    index = ["1"] * nframes
    with open(output, "w") as f:
        f.write(" ".join(index))


if __name__ == "__main__":
    reslice_dti()
