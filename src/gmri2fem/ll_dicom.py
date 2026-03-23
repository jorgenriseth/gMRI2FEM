import shutil
import subprocess
import tempfile
from pathlib import Path

import click
import numpy as np
import pydicom
from pydicom.errors import InvalidDicomError
import simple_mri as sm


def read_dicom_trigger_times(dicomfile):
    dcm = pydicom.dcmread(dicomfile)
    all_frame_times = [
        f.CardiacSynchronizationSequence[0].NominalCardiacTriggerDelayTime
        for f in dcm.PerFrameFunctionalGroupsSequence
    ]
    return np.unique(all_frame_times)

def read_legacy_dicom_trigger_times(dicomfile):
    dicomfile = Path(dicomfile)
    unique_trigger_times = set()
    possible_files = [p for p in dicomfile.parent.iterdir() if p.is_file()]
    for f in sorted(possible_files):
        try:
            ds = pydicom.dcmread(f, stop_before_pixels=True)
            unique_trigger_times.add(ds.TriggerTime)
        except (InvalidDicomError, AttributeError) as e:
            continue
    trigger_times = sorted(unique_trigger_times)
    if len(trigger_times) == 0:
        raise ValueError(f"Couldn't find any DICOM file with TriggerTime in {dicomfile.parent}")
    return trigger_times


def dcm2nii_looklocker(dicomfile, outpath):
    outdir, form = outpath.parent, outpath.stem
    outdir.mkdir(exist_ok=True, parents=True)
    try:
        times = read_dicom_trigger_times(dicomfile)
    except AttributeError as e:
        times = read_legacy_dicom_trigger_times(dicomfile)

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
        np.savetxt(f"{outdir}/{form}" + "_trigger_times.txt", times)


@click.command()
@click.option("--dicomfile", type=Path, required=True)
@click.option("--outpath", type=Path, required=True)
def dcm2nii_looklocker_cli(*args, **kwargs):
    dcm2nii_looklocker(*args, **kwargs)
