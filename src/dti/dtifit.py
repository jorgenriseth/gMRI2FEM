import subprocess
import tempfile
from pathlib import Path
from re import subn
from typing import Optional

from dti.utils import create_mask, with_suffix


def dtifit(
    dti_eddy_corrected: Path,
    bvals: Path,
    output: Path,
    tmppath: Optional[Path] = None,
):
    if tmppath is None:
        tmpdir = tempfile.TemporaryDirectory()
        tmppath = Path(tmpdir.name)

    eddy_b0 = tmppath / "eddy_corrected_b0.nii.gz"
    b0_cmd = f"fslroi {dti_eddy_corrected} {eddy_b0} 150 10"
    subprocess.run(b0_cmd, shell=True, check=True)

    eddy_b0_mean = tmppath / "eddy_corrected_b0_mean.nii.gz"
    mean_cmd = f"fslmaths {eddy_b0} -Tmean {eddy_b0_mean}"
    subprocess.run(mean_cmd, shell=True, check=True)

    mask = tmppath / "eddy_mask.nii.gz"
    create_mask(
        eddy_b0_mean,
        mask.parent / mask.name.replace("_mask.nii.gz", ""),
        threshold=0.4,
    )

    bvecs = with_suffix(dti_eddy_corrected, ".eddy_rotated_bvecs")
    fit_cmd = (
        "dtifit"
        + f" -k {dti_eddy_corrected}"
        + f" -o {with_suffix(output, '')}"
        + f" -m {mask}"
        + f" -r {bvecs}"
        + f" -b {bvals}"
        + " --save_tensor"
    )
    subprocess.run(fit_cmd, shell=True, check=True)
