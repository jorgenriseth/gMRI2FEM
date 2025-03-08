import json
import subprocess
from pathlib import Path


def mri_number_of_frames(input: str | Path) -> int:
    return int(
        subprocess.check_output(
            f"mri_info --nframes {input} | grep -v -E 'INFO|unknown time'", shell=True
        )
    )


def readout_time(sidecar: Path) -> str:
    with open(sidecar, "r") as f:
        meta = json.load(f)
    return meta["EstimatedTotalReadoutTime"]


def with_suffix(p: Path, newsuffix: str) -> Path:
    return p.parent / f"{p.name.split('.')[0]}{newsuffix}"
