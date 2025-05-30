import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import click
import tqdm

from dti.utils import mri_number_of_frames


# TODO: Incorporate to allow multiple files resliced simultaneously
def moving_pair(args):
    assert len(args) == 2, f"Got {len(args)} arguments, should only have two"
    return f" -rm {args[0]} {args[1]} "


@click.command()
@click.option("--inpath", type=Path, required=True)
@click.option("--fixed", type=Path, required=True)
@click.option("--outpath", type=Path, required=True)
@click.option("--transform", type=Path)
@click.option("--threads", type=int, default=1)
@click.option("--tmppath", type=Path)
@click.option("--greedyargs", type=str)
def reslice4d(**kwargs):
    reslice_4d(**kwargs)


def reslice_4d(
    inpath: Path,
    fixed: Path,
    outpath: Path,
    transform: Optional[Path] = None,
    threads: int = 1,
    tmppath: Optional[Path] = None,
    greedyargs: Optional[str] = None,
) -> Path:
    if transform is None:
        transform = Path("")
    if tmppath is None:
        tmpfile = tempfile.TemporaryDirectory(prefix=Path(outpath).stem)
        tmppath = Path(tmpfile.name)
    if greedyargs is None:
        greedyargs = ""

    nframes = mri_number_of_frames(inpath)

    for i in range(nframes):
        tmp_split = tmppath / f"slice{i}.nii.gz"
        tmp_reslice = tmppath / f"reslice{i}.nii.gz"
        subprocess.run(
            f"fslroi {inpath} {tmp_split} {i} 1", shell=True
        ).check_returncode()
        subprocess.run(
            f"greedy -d 3 -rf {fixed} {greedyargs} -rm {tmp_split} {tmp_reslice} -r {transform} -threads {threads} ",
            shell=True,
        ).check_returncode()
    components = [str(tmppath / f"reslice{i}.nii.gz") for i in range(nframes)]
    try: 
        subprocess.run(
            f"fslmerge -t {outpath} {' '.join(components)}", shell=True
        ).check_returncode()
    except subprocess.CalledProcessError:
        # Potential error due to slow I/O, wait a bit and retry.
        import time
        time.sleep(10)
        subprocess.run(
            f"fslmerge -t {outpath} {' '.join(components)}", shell=True
        ).check_returncode()
    return outpath
