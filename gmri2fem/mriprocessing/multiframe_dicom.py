import json
import logging
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional
from uuid import uuid4

import pydicom

from gmri2fem.filters import is_datedir, is_niifile
from gmri2fem.settings import DICOMSettings, patient_data_settings

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def patient_dicom2nii(
    dicom_patientdir: Path,
    t1dir: Path,
    t2dir: Path,
    looklockerdir: Path,
    patterns: dict[str, str],
    dtidir: Optional[Path] = None,
):
    t2_idx = 0
    for idx, studydir in enumerate(study_iterator(dicom_patientdir)):
        dicom2nii(studydir, t1dir, "T1", patterns["T1"])
        dicom2nii(studydir, looklockerdir, "LookLocker", patterns["LookLocker"])
        if idx == t2_idx:
            found = dicom2nii(studydir, t2dir, "T2", patterns["T2"])
            if not found:
                t2_idx += 1
        if dtidir is not None:
            raise NotImplementedError("Need to find correct DTI sequence.")


def dicom2nii(dicom_studydir: Path, outputdir: Path, label: str, pattern: str):
    outputdir.mkdir(exist_ok=True, parents=True)
    tempdir = Path(f"/tmp/gmri2fem_convert-{dicom_studydir.stem}-{uuid4()}")
    tempdir.mkdir()

    all = True

    try:
        path = find_dicom_sequence_dir(dicom_studydir, pattern)
        cmd = (
            "dcm2niix"
            + f" -f {label}"
            + f" -o '{tempdir}'"
            + f" '{path}/DICOM'"
            + f" 1>> '{tempdir}/log.txt'"
        )
        logger.info(f"Running conversion command: {cmd}")
        subprocess.run(cmd, shell=True).check_returncode()
    except ValueError as e:
        if "No DICOM" in str(e):
            logger.warn(str(e))
        all = False

    for nii_file in filter(is_niifile, tempdir.iterdir()):
        datedir = dicom_studydir.parent
        with open(nii_file.with_suffix(".json")) as f:
            info = json.load(f)
        timestamp = datetime.strptime(info["AcquisitionTime"], "%H:%M:%S.%f").strftime(
            "%H%M%S"
        )
        outputname = f"{datedir.stem.replace('_', '')}_{timestamp}"
        if "LookLocker" in nii_file.stem:
            outputname += f"_t{info['TriggerDelayTime']}"

        nii_file.rename(outputdir / f"{outputname}.nii")

    return all


def find_dicom_sequence_dir(dicom_studydir: Path, pattern: str) -> Path:
    paths = [
        x for x in sorted(dicom_studydir.glob("DICOM_*")) if re.match(pattern, x.name)
    ]
    if len(paths) == 0:
        raise ValueError(
            f"No DICOM directory found in {dicom_studydir} matching '{pattern}'."
        )
    if len(paths) > 1:
        warn = f"Multiple DICOM directories found in {dicom_studydir} matching '{pattern}': {paths}"
        logger.warning(warn)
    return paths[-1]


def find_timestamp_from_dicom(dicom_file: Path) -> Optional[str]:
    with pydicom.dcmread(dicom_file) as f:
        timestamp = f.StudyTime
    return timestamp


def study_iterator(dicom_patientdir: Path) -> list[Path]:
    return [
        study
        for date in sorted(filter(is_datedir, dicom_patientdir.iterdir()))
        for study in sorted(date.iterdir())
    ]


if __name__ == "__main__":
    from argparse import ArgumentParser

    patient_parser = ArgumentParser()
    patient_parser.add_argument(
        "patientid", type=str, help="Patient ID on format sub-XX"
    )
    patient_args, remaining_args = patient_parser.parse_known_args()
    paths = patient_data_settings(patientid=patient_args.patientid)
    settings = DICOMSettings(paths=paths)

    """Possible to override defaults."""
    parser = ArgumentParser()
    parser.add_argument("--dicomdir", type=str, default=paths.dicompath)
    parser.add_argument("--t1_outputdir", type=str, default=paths.t1raw)
    parser.add_argument("--t2_outputdir", type=str, default=paths.t2)
    parser.add_argument("--ll_outputdir", type=str, default=paths.looklocker)
    args = parser.parse_args(remaining_args)

    patient_dicom2nii(
        dicom_patientdir=Path(args.dicomdir),
        t1dir=args.t1_outputdir,
        t2dir=args.t2_outputdir,
        looklockerdir=args.ll_outputdir,
        patterns=settings.patterns,
    )

    for path in Path("/tmp").iterdir():
        if "gmri2fem_convert" in path.stem:
            shutil.rmtree(path)
