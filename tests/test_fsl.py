import os
from pathlib import Path
import shutil
import subprocess

import nibabel as nib
import numpy as np
import pytest


@pytest.mark.parametrize("command", ["eddy", "topup", "dtifit", "bet", "fslroi", "fslmerge", "fslmaths"])
def test_fsl_executable(command):
    assert shutil.which(command), f"Missing FSL executable: {command}"


def test_topup_configuration():
    assert (Path(os.environ["FSLDIR"]) / "etc/flirtsch/b02b0.cnf").is_file()


def test_fsl_utilities(tmp_path):
    source = tmp_path / "input.nii.gz"
    roi = tmp_path / "roi.nii.gz"
    merged = tmp_path / "merged.nii.gz"
    mean = tmp_path / "mean.nii.gz"
    data = np.ones((3, 3, 3, 2), dtype=np.float32)
    nib.save(nib.Nifti1Image(data, np.eye(4)), source)
    subprocess.run(["fslroi", str(source), str(roi), "0", "1"], check=True)
    subprocess.run(["fslmerge", "-t", str(merged), str(roi), str(roi)], check=True)
    subprocess.run(["fslmaths", str(merged), "-Tmean", str(mean)], check=True)
    assert nib.load(merged).shape == data.shape
    assert np.allclose(nib.load(mean).get_fdata(), 1)
