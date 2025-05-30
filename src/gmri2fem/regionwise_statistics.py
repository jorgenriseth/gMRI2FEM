from pathlib import Path

import numpy as np
import pandas as pd
import tqdm
from simple_mri import assert_same_space, load_mri

from gmri2fem.segmentation_groups import default_segmentation_groups
from gmri2fem.segment_tools import read_lut, find_label_description

REGION_STAT_FUNCTIONS = {
    "mean": np.mean,
    "median": np.median,
    "std": np.std,
    "PC1": lambda x: np.quantile(x, 0.01),
    "PC5": lambda x: np.quantile(x, 0.05),
    "PC95": lambda x: np.quantile(x, 0.95),
    "PC99": lambda x: np.quantile(x, 0.99),
}


def find_timestamp(
    timetable_path: Path, timestamp_sequence: str, subject: str, session: str
) -> float:
    try:
        timetable = pd.read_csv(timetable_path, sep="\t")
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
    timestamps_path: Path,
) -> pd.DataFrame:
    seg_mri = load_mri(seg_path, dtype=np.int16)
    data_mri = load_mri(mri_path, dtype=np.single)
    assert_same_space(seg_mri, data_mri)
    seg, data = seg_mri.data, data_mri.data

    seg_labels = np.unique(seg[seg != 0])
    fs_lut = read_lut(None)
    seg_descriptions = [
        fs_lut["description"][fs_lut["label"] == label].item() for label in seg_labels
    ]
    regions = {
        **{region: [int(label)] for region, label in zip(seg_descriptions, seg_labels)},
        **default_segmentation_groups(),
    }
    dframe = compute_region_statistics(data, seg, regions)

    timestamp = find_timestamp(timestamps_path, timestamp_sequence, subject, session)
    timestamp = max(0, timestamp)

    new_cols = ["subject", "session", "sequence", "timestamp"] + list(dframe.columns)
    newframe = dframe.assign(
        subject=subject, session=session, sequence=sequence, timestamp=timestamp
    )
    return newframe.loc[:, new_cols]


def compute_region_statistics(
    volume: np.ndarray, seg_vol: np.ndarray, regions: dict[str, list[int]]
) -> pd.DataFrame:
    records = []
    finite_mask = np.isfinite(volume)
    for description, labels in tqdm.tqdm(regions.items()):
        region_mask = np.isin(seg_vol, labels) * finite_mask
        region_data = volume[region_mask]

        voxelcount = region_mask.sum()
        if voxelcount == 0:
            continue

        group_regions = {
            **{
                "FS_LUT-labels": ",".join([str(x) for x in labels]),
                "FS_LUT-region": description,
                "FS_LUT-voxelcount": voxelcount,
                "region_total": np.sum(region_data),
            },
            **{
                f"{stat}": func(region_data) if voxelcount > 0 else np.nan
                for stat, func in REGION_STAT_FUNCTIONS.items()
            },
        }
        records.append(group_regions)
    return pd.DataFrame.from_records(records)

def voxel_count_to_ml_scale(affine: np.ndarray):
    return 1e-3 * np.linalg.det(affine[:3, :3])

def segstats(seg: np.ndarray, lut: pd.DataFrame, volscale: float):
    labels = np.unique(seg[seg > 0])
    seg_table = pd.DataFrame.from_records([
        {
            "label": label,
            "description": find_label_description(label, lut),
            "voxelcount": (seg == label).sum(),
            "volume (mL)": volscale * (seg == label).sum() 
        } for label in labels
    ])
    total = {
        "label": set(labels),
        "description": "all-regions",
        "voxelcount": (seg != 0).sum(),
        "volume (mL)": volscale * (seg != 0).sum() 
    }
    seg_table = pd.concat([
            seg_table,
            pd.DataFrame.from_records([total]),
        ],
        ignore_index=True,
    )
    return seg_table


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--subject", type=str, required=True)
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
        subject=args.subject,
        session=args.subject_session,
        sequence=args.sequence,
        timestamp_sequence=args.timestamp_sequence,
        mri_path=args.data,
        seg_path=args.seg,
        timestamps_path=args.timestamps,
    )
    dframe.to_csv(Path(args.output), index=False, sep=";")
