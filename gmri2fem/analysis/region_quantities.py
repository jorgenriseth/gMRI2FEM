import json
from pathlib import Path

import nibabel.nifti1 as nifti1
import numpy as np
import pandas as pd
from tqdm import tqdm

from gmri2fem.analysis.seg_groups import default_segmentation_groups

# Intervals used to group the T1mapvarying timestamps of different subjects
# together. Intervals are given in hours, and data gets the label 't{i}'
# i.e. t=4.5h -> 't1', t=23.2h -> 't2', and so on.
SESSION_INTERVALS = [(0, 0), (2, 6), (20, 28), (45, 52), (68, 80)]


def timestamp_session_idx(t: float) -> int:
    """Get the index of the session corresponding to the given timestamp."""
    for idx, interval in enumerate(SESSION_INTERVALS):
        if interval[0] <= t / 3600 <= interval[1]:
            return idx + 1
    raise ValueError(f"Timestamp {t}s={t/3600:.2f}h is not in any session intervals")


def find_timestamp(
    timetable_path: Path, timestamp_sequence: str, subject: str, session: str
) -> float:
    try:
        timetable = pd.read_csv(timetable_path)
    except pd.errors.EmptyDataError:
        raise RuntimeError(f"Timetable-file {timetable_path} is empty.")

    try:
        timestamp = timetable.loc[
            (timetable["sequence_label"] == timestamp_sequence)
            & (timetable["subject"] == subject)
            & (timetable["session"] == session)
        ]["acquisition_relative_injection"]
    except ValueError as e:
        print(timetable)
        print(timestamp_sequence, subject, session)
        raise e
    return timestamp.item()


def create_dataframe(
    subject: str,
    session: str,
    sequence: str,
    timestamp_sequence: str,
    mri_path: Path,
    seg_path: Path,
    lut_path: Path,
    timestamps_path: Path,
) -> pd.DataFrame:
    with open(lut_path, "r") as f:
        fs_lut = json.load(f)

    seg = nifti1.load(seg_path).get_fdata(dtype=np.single).astype(np.int32)

    # Extract all freesurfer-lookuptable values and region names present in segmentation
    seg_labels = np.unique(seg[seg != 0])
    seg_regions = [fs_lut[str(label)] for label in seg_labels]

    regions = {
        **{region: [label] for region, label in zip(seg_regions, seg_labels)},
        **default_segmentation_groups(),
    }

    data = nifti1.load(mri_path).get_fdata(dtype=np.single)
    timetable = pd.read_csv(timestamps_path)
    timestamp = find_timestamp(timestamps_path, timestamp_sequence, subject, session)

    regions_stat_functions = {
        "mean": np.mean,
        "median": np.median,
        "std": np.std,
        "PC1": lambda x: np.quantile(x, 0.01),
        "PC5": lambda x: np.quantile(x, 0.05),
        "PC95": lambda x: np.quantile(x, 0.95),
        "PC99": lambda x: np.quantile(x, 0.99),
    }

    timestamp = max(0, timestamp)
    # In case subject missed a session, the session_idx might not represent the
    # corresponding session for the rest of the subjects. Therefore, identify
    # session-idx based on timestamp and interval instead.
    session_idx = timestamp_session_idx(timestamp)

    records = []
    for name, labels in tqdm(regions.items()):
        regions_data = data[np.isin(seg, labels)]
        regions_data = regions_data[~np.isnan(regions_data)]  # Eliminate nans

        if regions_data.size == 0:
            continue

        group_regions = {
            **{
                "FS_LUT-labels": ",".join([str(x) for x in labels]),
                "FS_LUT-region": name,
                "FS_LUT-voxelcount": seg[np.isin(seg, labels)].size,
                "subject": subject,
                "session": f"ses-{session_idx:02d}",
                "sequence": sequence,
                "time": timestamp,
                "region_total": np.sum(regions_data),
            },
            **{
                f"{stat}": func(regions_data)
                for stat, func in regions_stat_functions.items()
            },
        }
        records.append(group_regions)

    return pd.DataFrame.from_records(records)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--subjectid", type=str, required=True)
    parser.add_argument("--subject_session", type=str, required=True)
    parser.add_argument("--sequence", type=str, required=True)
    parser.add_argument("--timestamp_sequence", type=str)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--seg", type=Path, required=True)
    parser.add_argument("--timestamps", type=Path, required=True)
    parser.add_argument("--lutfile", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    if args.timestamp_sequence is None:
        args.timestamp_sequence = args.sequence
    dframe = create_dataframe(
        subject=args.subjectid,
        session=args.subject_session,
        sequence=args.sequence,
        timestamp_sequence=args.timestamp_sequence,
        mri_path=args.data,
        seg_path=args.seg,
        lut_path=args.lutfile,
        timestamps_path=args.timestamps,
    )
    dframe.to_csv(Path(args.output), index=False)
