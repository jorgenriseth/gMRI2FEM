import subprocess
from pathlib import Path


def mri_number_of_frames(input: str | Path) -> int:
    return int(
        subprocess.check_output(
            f"mri_info --nframes {input} | grep -v -E 'INFO|unknown time'", shell=True
        )
    )


def with_suffix(p: Path, newsuffix: str) -> Path:
    return p.parent / f"{p.name.split('.')[0]}{newsuffix}"


def path_stem(p: Path) -> str:
    """Returns path stem, keeping only what is before the first dot."""
    return f"{p.name.split('.')[0]}"


def create_mask(input: Path, output: Path, threshold: float):
    mask_cmd = f"bet {input} {output.parent / path_stem(output)} -m -f {threshold} -n"
    subprocess.run(mask_cmd, shell=True, check=True)
