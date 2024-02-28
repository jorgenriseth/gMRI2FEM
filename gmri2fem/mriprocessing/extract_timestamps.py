"""
Script to extract the time in seconds since tracer injection in for each subject
"""
import datetime
import json
from pathlib import Path

import numpy as np

from gmri2fem.filters import is_T1_mgz


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


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--subject", type=str, required=True)
    parser.add_argument("--injection_times", type=Path, required=True)
    parser.add_argument("--t1dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(list(Path(args.t1dir).iterdir()))

    times = timestamps(
        args.subject,
        Path(args.injection_times),
        sorted(filter(is_T1_mgz, Path(args.t1dir).iterdir())),
    )
    print(times)
    np.savetxt(args.output, times)
