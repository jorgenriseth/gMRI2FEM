import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from dti.utils import create_mask, mri_number_of_frames
from gmri2fem.utils import with_suffix


def eddy_correct(
    dti,
    topup_b0_mean,
    acq_params,
    output: Path,
    multiband_factor: int = 1,
    nthreads: int = 1,
    tmppath: Optional[Path] = None,
    verbose: bool = False,
):
    if tmppath is None:
        tmpdir = tempfile.TemporaryDirectory()
        tmppath = Path(tmpdir.name)

    index_file = tmppath / "eddy_index.txt"
    create_eddy_index_file(dti, index_file)

    mask = tmppath / "topup_mask.nii.gz"
    create_mask(
        topup_b0_mean, Path(str(mask).replace("_mask.nii.gz", "")), threshold=0.2
    )

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
        + f" --nthr={nthreads}"
        + " --verbose" * verbose
    )
    subprocess.run(eddy_cmd, shell=True, check=True)


def create_eddy_index_file(input: Path, output: Path):
    index = ["1"] * mri_number_of_frames(input)
    with open(output, "w") as f:
        f.write(" ".join(index))
