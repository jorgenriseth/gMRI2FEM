import subprocess
import nibabel
from nibabel.spatialimages import SpatialImage
from pathlib import Path


def load_nibabel_spatial(input: str | Path) -> SpatialImage:
    im = nibabel.load(input)
    assert isinstance(im, SpatialImage), f"File {input} is not a nibabel SpatialImage."
    return im


def mri_number_of_frames(input: str | Path) -> int:
    im = load_nibabel_spatial(input)
    ndim = im.ndim
    if ndim == 3:
        return 1
    elif ndim == 4:
        return im.shape[3]
    raise RuntimeError(f"Invalid dimension {ndim}, of mri {input}")


def path_stem(p: Path) -> str:
    """Returns path stem, keeping only what is before the first dot."""
    return f"{p.name.split('.')[0]}"


def create_mask(input: Path, output: Path, threshold: float):
    mask_cmd = f"bet {input} {output.parent / path_stem(output)} -m -f {threshold} -n"
    subprocess.run(mask_cmd, shell=True, check=True)
