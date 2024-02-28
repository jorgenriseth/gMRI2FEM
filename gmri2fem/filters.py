import re
from pathlib import Path


def is_not_template(p: Path) -> bool:
    return p.suffix == ".mgz" and p.stem != "template"


def is_datetime_dir(path: Path) -> bool:
    match = re.match("[0-9]{8}_[0-9]{6}$", path.name)
    return path.is_dir() and (match is not None)


def is_datedir(path: Path) -> bool:
    match = re.match("[0-9]{4}_[0-9]{2}_[0-9]{2}$", path.name)
    return path.is_dir() and (match is not None)


def is_niifile(path: Path) -> bool:
    return path.suffix == ".nii"


def is_T1_mgz(p: Path) -> bool:
    return re.match("[0-9]{8}_[0-9]{6}.mgz$", p.name) is not None


def is_datetime_nii(p: Path):
    return re.match("[0-9]{8}_[0-9]{6}.nii$", p.name) is not None
