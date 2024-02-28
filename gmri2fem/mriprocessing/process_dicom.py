"""Program to run through the DICOM-file-structure that is typically given to 
us from the hospital, and extract the files of interesest.
"""
import csv
import datetime
import json
import re
import shutil
import subprocess
from pathlib import Path

import pydicom
from loguru import logger


def generate_metadata(
    participants_private: Path, participants_tsv: Path, sourcedata: Path
):
    with open(participants_private, "r") as f:
        private_data = json.load(f)
    with open("participants.json", "r") as f:
        header = ["participant_id", *json.load(f).keys()]

    with open(participants_tsv, "w") as f:
        w = csv.DictWriter(f, header, delimiter="\t")
        w.writeheader()
        for subject in private_data:
            record = participant_record(subject, sourcedata / subject["id"])
            w.writerow(record)
    return private_data


def participant_record(
    subject: dict[str, str], subject_dicomdir: Path
) -> dict[str, str]:
    participant_record = {
        "participant_id": subject["id"],
        "age": subject["age"],
        "sex": subject["sex"],
    }
    session_dicomdirs = subject_session_dirs(subject_dicomdir)
    assert len(session_dicomdirs) > 0
    for idx, session in enumerate(session_dicomdirs):
        session_time = session_datetime(session)
        injection_time = datetime.datetime.strptime(
            subject["injection_time"], "%Y%m%d_%H%M%S"
        )
        dt = session_time.date() - injection_time.date()
        participant_record[f"ses-{idx+1:02d}-day"] = str(dt.days)
        participant_record[f"ses-{idx+1:02d}-tod"] = session_time.strftime("%H:%M:%S")
    return participant_record


def session_datetime(session_dicomdir: Path) -> datetime.datetime:
    anysequence = next(session_dicomdir.rglob("IM_0002"))
    with pydicom.dcmread(anysequence) as ds:
        session_time = ds.StudyTime
        session_date = ds.StudyDate  # Probably need to switch
    return datetime.datetime.combine(
        datetime.datetime.strptime(session_date, "%Y%m%d").date(),
        datetime.datetime.strptime(session_time, "%H%M%S").time(),
    )


def process_subject_dicom(
    subject_id: str,
    subject_dicomdir: Path,
    subject_dicomdir_out: Path,
    subject_niftidir: Path,
    patterns: dict[str, str],
    ignore_existing: bool,
):
    sessions = subject_session_dirs(subject_dicomdir)
    for label, pattern in patterns.items():
        sequence_dirs = sequence_directories(sessions, label, pattern)
        for idx, sequencedir in enumerate(sequence_dirs):
            session = f"{idx+1:02d}"
            subdir = bids_subdir_pattern(session, label)
            copy_dicom_sequence(
                sequencedir,
                subject_dicomdir_out / f"{subdir}/{label}",
                ignore_existing,
            )
            dicom2nii(
                sequencedir,
                subject_niftidir / subdir,
                f"sub-{subject_id}_ses-{session}_{label}",
            )


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
        shutil.copytree(sequencedir_raw, sequencedir_out)
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
        cmd = f"dcm2niix -w 2 -z y -f {label} -o '{outputdir}' '{sequencedir}/DICOM' >> log_dcm2niix.txt "
        logger.info(f"Running conversion command: {cmd}")
        subprocess.run(cmd, shell=True).check_returncode()
    except ValueError as e:
        if "No DICOM" in str(e):
            logger.warning(str(e))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--sourcedata", type=Path, default="sourcedata/")
    parser.add_argument("--curated", type=Path, default="sourcedata-curated/")
    parser.add_argument("--ignore_existing", action="store_true")
    args = parser.parse_args()

    patterns = {
        "T1w": r"(DICOM_\d+_\d+[_ ])(.*(T1_3D).*)",
        "T2w": r"(DICOM_\d+_\d+[_ ])(.*(TE565|T2W).*)",
        "LookLocker": r"(DICOM_\d+_3[_ ])(.*(LookLocker|2beatpause).*)",
        "Mixed": r"(DICOM_\d+_\d+[_ ])(.*(Mixed).*)",
        "FLAIR": r"(DICOM_\d+_\d+[_ ])(.*(FLAIR 3D).*)",
    }
    subjects = generate_metadata(
        participants_private=Path("participants-private.json"),
        participants_tsv=Path("participants.tsv"),
        sourcedata=args.sourcedata,
    )
    for subject in subjects:
        subject_raw_dicomdir = args.sourcedata / str(subject["id"])
        subject_structured_dicomdir = args.curated / str(subject["id"])
        subject_niftidir = Path(str(subject["id"]))
        try:
            process_subject_dicom(
                subject["id"],
                subject_raw_dicomdir,
                subject_structured_dicomdir,
                subject_niftidir,
                patterns,
                args.ignore_existing,
            )
        except ValueError as e:
            if "in any of the sessions" in str(e):
                logger.warning(e)
            else:
                raise e
