"""
Script to extract the time in seconds since tracer injection in for each subject
"""
import datetime
import json
import pandas as pd
from pathlib import Path

import numpy as np


def image_timestamp(p: Path) -> datetime.datetime:
    return datetime.datetime.strptime(p.stem, "%Y%m%d_%H%M%S")


def injection_timestamp(injection_time_file: Path, subjectid: str) -> datetime.datetime:
    with open(injection_time_file, "r") as f:
        time_string = json.load(f)[subjectid]
    return datetime.datetime.strptime(time_string.strip(), "%Y%m%d_%H%M%S")


def timestamps(subject, timestamp_file, t1_files):
    t0 = injection_timestamp(timestamp_file, subject)
    times = np.zeros(len(t1_files))
    print(t0, subject, timestamp_file, t1_files)
    for idx, cpath in enumerate(t1_files):
        t_idx = image_timestamp(cpath)
        times[idx] = max(0, int(np.rint((t_idx - t0).total_seconds())))
    return times


def read_timetable(
    timetable_path: Path, subjectid: str, sequence_label: str
) -> np.ndarray:
    dframe = pd.read_csv(
        timetable_path,
    )
    subject_sequence_entries = (dframe["subject"] == subjectid) & (
        dframe["sequence_label"] == sequence_label
    )
    acq_times = dframe[subject_sequence_entries]["acquisition_relative_injection"]
    return np.array(acq_times)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--timetable", type=Path, required=True)
    parser.add_argument("--subject", type=str, required=True)
    parser.add_argument("--sequence_label", type=str, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    times = read_timetable(args.timetable, args.subject, args.sequence_label)
    assert len(times) > 0
    np.savetxt(args.output, np.maximum(0.0, times))
