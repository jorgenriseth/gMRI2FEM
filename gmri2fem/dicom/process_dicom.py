"""Program to run through the DICOM-file-structure that is typically given to 
us from the hospital, and extract the files of interesest.
"""
import json
import re
import shutil
import subprocess
import datetime
from pathlib import Path

import nibabel
import pandas as pd
import pydicom
from loguru import logger

from gmri2fem.dicom.extract_mixed_sequence import estimate_T1_mixed, extract_volume


MIXED_VOLUME_LABELS = [
    "IR-modulus",
    "IR-real",
    "IR-phase_corrected",
    "SE-modulus",
    "SE-real",
    "T1-scanner",
]


def process_subject_dicom(
    subject: dict[str, str],
    sourcedata: Path,
    sourcedata_curated: Path,
    derived: Path,
    patterns: dict[str, str],
    ignore_existing: bool,
):
    subject_id = subject["subject-id"]
    subject_dicomdir = sourcedata / subject["subject-id"]
    subject_dicomdir_out = sourcedata_curated / subject["subject-id"]
    subject_niftidir = Path(subject["subject-id"])
    subject_derived = derived / subject["subject-id"]
    sessions = subject_session_dirs(subject_dicomdir)
    for label, pattern in patterns.items():
        sequence_dirs = sequence_directories(sessions, label, pattern)
        for idx, sequencedir in enumerate(sequence_dirs):
            session = f"{idx+1:02d}"
            subdir = bids_subdir_pattern(session, label)
            sequencedir_out = subject_dicomdir_out / f"{subdir}/{label}"
            try:
                shutil.copytree(sequencedir, sequencedir_out)
            except FileExistsError as e:
                if ignore_existing:
                    logger.info(f"{sequencedir_out} already exists. Continuing.")
                    continue
                raise e
            dicom2nii(
                sequencedir,
                subject_niftidir / subdir,
                f"{subject_id}_ses-{session}_{label}",
            )
            if "Mixed" in label:
                mixed_paths = sorted(sequencedir.glob("DICOM/IM_*"))
                assert (
                    len(mixed_paths) == 1
                ), f"Couldn't find single IM-file in {sequencedir}"
                mixed_dir_out = subject_derived / f"{subdir}"
                mixed_dir_out.mkdir(exist_ok=True, parents=True)
                filename = f"{subject_id}_ses-{session}_T1map_mixed.nii.gz"
                T1map = estimate_T1_mixed(mixed_paths[0], T1_lo=200, T1_hi=5000)
                nibabel.nifti1.save(T1map, mixed_dir_out / filename)


def subject_session_dirs(subject_dicomdir: Path) -> list[Path]:
    return [
        session_dir
        for date_dir in sorted(subject_dicomdir.glob("202[34]*"))
        for session_dir in sorted(date_dir.iterdir())
    ]


def find_sequence_in_session(session_dicomdir: Path, pattern: str) -> Path:
    paths = [
        x for x in sorted(session_dicomdir.glob("DICOM_*")) if re.match(pattern, x.name)
    ]
    if len(paths) == 0:
        raise ValueError(
            f"No DICOM directory found in {session_dicomdir} matching '{pattern}'."
        )
    if len(paths) > 1:
        warn = f"Multiple DICOM directories found in {session_dicomdir} matching '{pattern}':\n"
        warn += "\n".join([str(p) for p in paths]) + "\n"
        warn += "Returning final instance: " + str(paths[-1])
        logger.warning(warn)
    return paths[-1]


def sequence_directories(sessions: list[Path], label: str, pattern: str) -> list[Path]:
    sequence_dirs = []
    for sessiondir in sessions:
        try:
            sequence_dirs.append(find_sequence_in_session(sessiondir, pattern))
        except ValueError:
            pass
    validate_sequence_dirs(sessions, sequence_dirs, label, pattern)
    return sequence_dirs


def validate_sequence_dirs(
    sessions: list[Path], sequence_dirs: list[Path], label: str, pattern: str
):
    if label in ["T2w", "FLAIR"]:
        if len(sequence_dirs) == 0:
            e = f"Couldn't find pattern {pattern} in any of the sessions:\n"
            e += "\n".join([str(p) for p in sessions])
            raise ValueError(e)
        elif len(sequence_dirs) > 1:
            e = f"Found pattern {pattern} in multiple sessions:\n"
            e += "\n".join([str(p) for p in sessions])
            raise ValueError(e)
    elif label in ["T1w", "LookLocker", "Mixed"]:
        assert len(sequence_dirs) == len(sessions)
    else:
        raise ValueError(f"Unrecognized label {label}, with pattern {pattern}")


def bids_subdir_pattern(session: str, label: str) -> str:
    if label in ["T1w", "LookLocker", "Mixed", "T2w", "FLAIR"]:
        return f"ses-{session}/anat/"
    elif label in ["DTI"]:
        return f"ses-{session}/dwi/"
    else:
        raise ValueError(f"No rule for placing sequence labeled {label}")


def copy_dicom_sequence(
    sequencedir_raw: Path, sequencedir_out: Path, ignore_existing=True
):
    try:
        logger.info(f"Successfully copied {sequencedir_raw} to {sequencedir_out}")
        return True
    except FileExistsError as e:
        if ignore_existing:
            pass
            logger.info(f"{sequencedir_out} already exists. Continuing.")
        else:
            raise e


def dicom2nii(sequencedir: Path, outputdir: Path, label: str):
    outputdir.mkdir(exist_ok=True, parents=True)
    try:
        cmd = f"dcm2niix -w 2 -z y -f {label} -o '{outputdir}' '{sequencedir}/DICOM' >> /tmp/dcm2niix.txt "
        logger.info(f"Running conversion command: {cmd}")
        subprocess.run(cmd, shell=True).check_returncode()
    except ValueError as e:
        if "No DICOM" in str(e):
            logger.warning(str(e))

    if "Mixed" in label:
        logger.info(f"Extracting Mixed sequence {sequencedir} > {outputdir / label}")
        dicom2nii_mixed(sequencedir, outputdir, label)


def subject_session_dates(subject_dicomdir: Path) -> dict[str, datetime.datetime]:
    session_dirs = subject_session_dirs(subject_dicomdir)
    return {
        f"ses-{idx+1:02d}": datetime.datetime.strptime(
            session_dir.parent.stem.replace("_", ""), "%Y%m%d"
        )
        for idx, session_dir in enumerate(session_dirs)
    }


def acquisition_time(sidecar: Path) -> datetime.time:
    with open(sidecar, "r") as f:
        info = json.load(f)
    acq_time = datetime.datetime.strptime(info["AcquisitionTime"], "%H:%M:%S.%f")
    return acq_time.time()


def sequence_acquisition_record(
    subject: dict[str, str], labels: list[str]
) -> list[dict[str, str | float]]:
    subject_dicoms = Path("sourcedata/") / subject["subject-id"]
    injection_time = datetime.datetime.strptime(
        subject["injection_time"], "%Y%m%d_%H%M%S"
    )
    session_dates = subject_session_dates(subject_dicoms)
    subject_records = []
    for session, date in session_dates.items():
        for label in list({x: None for x in labels}.keys()):
            subject_nifti = Path(f'{subject["subject-id"]}/{session}')
            sequence_paths = subject_nifti.rglob(f"*{label}*.json")
            try:
                acq_time = acquisition_time(next(sequence_paths))
            except StopIteration:
                continue
            timedelta = datetime.datetime.combine(date, acq_time) - injection_time
            record = {
                "subject": subject["subject-id"],
                "session": session,
                "SequenceLabel": label,
                "acquisition_relative_injection": timedelta.total_seconds(),
            }
            subject_records.append(record)
    return subject_records


def dicom2nii_mixed(sequencedir: Path, outputdir: Path, label: str):
    mixed_paths = sorted(sequencedir.glob("DICOM/IM_*"))
    assert len(mixed_paths) == 1
    info = pydicom.dcmread(mixed_paths[0])
    for volume_idx, volume_label in enumerate(MIXED_VOLUME_LABELS):
        nii = extract_volume(info, volume_idx)
        nibabel.nifti1.save(nii, outputdir / f"{label}_{volume_label}.nii.gz")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--participants", type=Path, required=True)
    parser.add_argument("--sourcedata", type=Path, default="sourcedata/")
    parser.add_argument("--curated", type=Path, default="sourcedata-curated/")
    parser.add_argument("--derivatives", type=Path, default="derivatives/")
    parser.add_argument("--ignore_existing", action="store_true")
    args = parser.parse_args()
    args.ignore_existing = True

    with open(args.participants, "r") as f:
        subjects = json.load(f)

    patterns = {
        "T1w": r"(DICOM_\d+_\d+[_ ])(.*(T1_3D).*)",
        "T2w": r"(DICOM_\d+_\d+[_ ])(.*(TE565|T2W).*)",
        "LookLocker": r"(DICOM_\d+_3[_ ])(.*(LookLocker|2beatpause).*)",
        "Mixed": r"(DICOM_\d+_\d+[_ ])(.*(Mixed).*)",
        "FLAIR": r"(DICOM_\d+_\d+[_ ])(.*(FLAIR 3D).*)",
    }
    for subject in subjects:
        try:
            process_subject_dicom(
                subject,
                args.sourcedata,
                args.curated,
                args.derived,
                patterns,
                args.ignore_existing,
            )
        except ValueError as e:
            if "in any of the sessions" in str(e):
                logger.warning(e)
            else:
                raise e

    records = sum(
        (
            sequence_acquisition_record(subject, list(patterns.keys()))
            for subject in subjects
        ),
        start=[],
    )
    dframe = pd.DataFrame.from_records(records)
    dframe.to_csv(args.derivatives / "timetable.csv")
