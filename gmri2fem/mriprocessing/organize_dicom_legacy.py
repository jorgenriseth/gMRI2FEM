#!/usr/bin/env/python
import json
import logging
import shutil
from datetime import datetime
from functools import partial
from itertools import chain, repeat
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Dict, Iterator, Tuple

import pydicom

logger = logging.getLogger(__name__)


def main(
    input_dir: Path,
    output_dir: Path,
    sequences: dict[str, str],
    add_metadata: bool = False,
) -> None:
    filemap = create_protocol_filemap(input_dir, output_dir, sequences)
    for src, dest in filemap.items():
        dest.parent.mkdir(exist_ok=True, parents=True)
        shutil.copy(src, dest)

    if add_metadata:
        add_sorted_metadata(output_dir)


def create_protocol_filemap(
    input_dir: Path, output_dir: Path, sequences: Dict[str, str], n_jobs=None
) -> Dict[Path, Path]:
    """Runs through the folder structure with the MRI data we receive from the
    hospital, and creates a dictionary mapping a src-path sorts them according
    to 'patient - study_datetime - protocol'. Runs in parallel for each of the studies.
    """
    pool = Pool(n_jobs)
    task = partial(
        create_study_filemap,
        output_dir=output_dir,
        sequences=sequences,
    )
    results = pool.map(task, study_iterator(input_dir))
    return dict(chain(*(x.items() for x in results)))


def create_study_filemap(
    study_dir: Path, output_dir: Path, sequences: Dict[str, str]
) -> Dict[Path, Path]:
    date_dir = study_dir.parent
    patient = date_dir.parent.stem
    study_data = study_dir / "DICOM" / "DICOM"
    date = datetime.strptime(date_dir.stem, "%Y_%m_%d").strftime("%Y%m%d")
    timestamp = find_timestamp(study_dir)
    study_target = output_dir / patient / f"{date}_{timestamp}"

    filemap = {}
    for imfile, offset in study_imfiles(study_data):
        with pydicom.dcmread(imfile) as f:
            file_protocol = f.ProtocolName
            if file_protocol in sequences:
                filemap[imfile] = (
                    study_target
                    / sequences[file_protocol]
                    / renumber_imfile(imfile, offset)
                )
    return filemap


def find_timestamp(study_dir: Path) -> str:
    study_data_path = study_dir / "DICOM" / "DICOM"
    first_imfile, _ = next(study_imfiles(study_data_path))
    with pydicom.dcmread(first_imfile) as f:
        timestamp = f.StudyTime
    return timestamp


def study_imfiles(study_data_path: Path) -> Iterator[Path]:
    return chain(
        study_imfiles_in_dicomdir(study_data_path),
        *(
            study_imfiles_in_dicom_subdir(subdir)
            for subdir in filter(is_image_subdirectory, study_data_path.iterdir())
        ),
    )


def study_imfiles_in_dicomdir(study_data_path: Path) -> Iterator[Tuple[Path, int]]:
    return zip(filter(is_imfile, study_data_path.iterdir()), repeat(0))


def study_imfiles_in_dicom_subdir(
    study_subdir_path: Path,
) -> Iterator[Tuple[Path, int]]:
    return zip(
        filter(is_imfile, study_subdir_path.iterdir()),
        repeat(int(study_subdir_path.stem)),
    )


def study_iterator(input_dir: Path) -> Iterator[Path]:
    return (
        x
        for patient in input_dir.iterdir()
        for date in patient.iterdir()
        for x in date.iterdir()
    )


def renumber_imfile(imfile: Path, offset: int) -> str:
    return f"IM_{int(imfile.stem.split('_')[1]) + offset * 2048:04d}"


def is_image_subdirectory(path: Path) -> bool:
    return path.is_dir() and path.stem.isdigit()


def is_imfile(path: Path) -> bool:
    """Checks if path is a dicom IM-file."""
    return path.is_file() and path.stem[:3] == "IM_"


def create_study_metadata(study_path: Path):
    study_metadata = {}
    for path in filter(lambda x: x.is_dir(), study_path.iterdir()):
        protocol_filelist = [p for p in path.iterdir() if is_imfile(p)]
        label_list = [int(x.stem.split("_")[1]) for x in protocol_filelist]
        with pydicom.dcmread(protocol_filelist[0]) as f:
            timestamp = f.SeriesTime

        study_metadata[path.stem] = {
            "min": min(label_list),
            "max": max(label_list),
            "num_files": len(protocol_filelist),
            "series_time": timestamp,
        }
    return study_metadata


def store_study_metadata(study_path: Path) -> None:
    filepath = study_path / "info.json"
    with open(filepath, "w") as f:
        json.dump(create_study_metadata(study_path), f, indent=4)


def add_sorted_metadata(output_dir: Path) -> None:
    reorganized_study_iterator = (
        study for patient in output_dir.iterdir() for study in patient.iterdir()
    )
    for study in reorganized_study_iterator:
        store_study_metadata(study)


if __name__ == "__main__":
    inputdir = Path("GRIP")
    outputdir = Path("GRIP_SORTED")
    sequences = {
        "WIP PDT1_3D 08mm": "T1",
        "WIP PDT1_3D 1mm": "T1",
        "WIP DelRec - WIP 2beatpause1mm 3000 HR 21": "Look-Locker",
        "WIP 07mmTE565 3D TSE": "T2",
        "WIP T2W 3D TSE TE565": "T2",
    }
    filemap = create_study_filemap(
        inputdir / "PAT_002" / "2023_02_13" / "Zi_121304", outputdir, sequences
    )
    for src, dest in filemap.items():
        dest.parent.mkdir(exist_ok=True, parents=True)
        shutil.copy(src, dest)
    add_sorted_metadata(outputdir)
