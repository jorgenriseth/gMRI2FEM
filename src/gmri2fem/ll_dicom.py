import shutil
import subprocess
import tempfile
from pathlib import Path

import click
import numpy as np
import pydicom
import simple_mri as sm


def read_dicom_trigger_times(dicomfile):
    dcm = pydicom.dcmread(dicomfile)
    all_frame_times = [
        f.CardiacSynchronizationSequence[0].NominalCardiacTriggerDelayTime
        for f in dcm.PerFrameFunctionalGroupsSequence
    ]
    return np.unique(all_frame_times)


def dcm2nii_looklocker(dicomfile, outpath):
    outdir, form = outpath.parent, outpath.stem
    outdir.mkdir(exist_ok=True, parents=True)
    times = read_dicom_trigger_times(dicomfile)
    np.savetxt(f"{outdir}/{form}" + "_trigger_times.txt", times)

    with tempfile.TemporaryDirectory(prefix=outpath.stem) as tmpdir:
        tmppath = Path(tmpdir)
        cmd = f"dcm2niix -f {form} -z y --ignore_trigger_times -o '{tmppath}' '{dicomfile}' > /tmp/dcm2niix.txt"
        subprocess.run(cmd, shell=True, check=True)
        shutil.copy(
            tmppath / f"{form}.json",
            outpath.with_suffix(".json"),
        )
        mri = sm.load_mri(tmppath / f"{form}.nii.gz", dtype=np.double)
        sm.save_mri(
            mri, outpath.with_suffix(".nii.gz"), dtype=np.single, intent_code=2001
        )


@click.command()
@click.option("--dicomfile", type=Path, required=True)
@click.option("--outpath", type=Path, required=True)
def dcm2nii_looklocker_cli(*args, **kwargs):
    dcm2nii_looklocker(*args, **kwargs)
