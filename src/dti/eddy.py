import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from dti.utils import mri_number_of_frames, with_suffix


def eddy_correct(
    dti,
    topup_b0_mean,
    acq_params,
    output: Path,
    multiband_factor: int = 1,
    tmppath: Optional[Path] = None,
):
    if tmppath is None:
        tmpdir = tempfile.TemporaryDirectory()
        tmppath = Path(tmpdir.name)

    index_file = tmppath / "eddy_index.txt"
    create_eddy_index_file(dti, index_file)

    mask = tmppath / "topup_mask.nii.gz"
    mask_cmd = (
        f"bet {topup_b0_mean} {str(mask).replace('_mask.nii.gz', '')} -m -f 0.2 -n"
    )
    subprocess.run(mask_cmd, shell=True, check=True)

    bvecs = with_suffix(dti, ".bvec")
    bvals = with_suffix(dti, ".bval")
    eddy_cmd = (
        "eddy diffusion"
        + f" --imain={dti}"
        + f" --mask={mask}"
        + f" --acqp={acq_params}"
        + f" --index={index_file}"
        + f" --bvecs={bvecs}"
        + f" --bvals={bvals}"
        + f" --topup={str(topup_b0_mean).replace('_b0_mean.nii.gz', '')}"
        + f" --out={output.parent / output.stem.split('.')[0]}"
        + f" --ol_type=both"
        + f" --mb={multiband_factor}"
    )
    subprocess.run(eddy_cmd, shell=True, check=True)


def create_eddy_index_file(input: Path, output: Path):
    index = ["1"] * mri_number_of_frames(input)
    with open(output, "w") as f:
        f.write(" ".join(index))
