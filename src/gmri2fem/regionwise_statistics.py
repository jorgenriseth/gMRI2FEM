import re
from pathlib import Path
from typing import Optional

import click
import numpy as np
import pandas as pd
import simple_mri as sm
import tqdm

import gmri2fem.segment_tools as segtools
from gmri2fem.segment_tools import find_label_description, read_lut
from gmri2fem.segmentation_groups import default_segmentation_groups
from gmri2fem.utils import find_timestamp, with_suffix, prepend_info, find_timestamps


@click.command("mristats")
@click.option("--segmentation", "-s", "seg_path", type=Path, required=True)
@click.option("--mri", "-m", "mri_paths", multiple=True, type=Path, required=True)
@click.option("--output", "-o", type=Path, required=True)
@click.option("--timetable", "-t", type=Path)
@click.option("--timelabel", "-l", "timetable_sequence", type=str)
def compute_mri_stats(
    seg_path: str | Path,
    mri_paths: tuple[str | Path],
    output: str | Path,
    timetable: Optional[str | Path],
    timetable_sequence: Optional[str | Path],
):
    if not Path(seg_path).exists():
        raise RuntimeError(f"Missing segmentation: {seg_path}")

    for path in mri_paths:
        if not Path(path).exists():
            raise RuntimeError(f"Missing: {path}")

    dataframes = [
        create_dataframe(
            Path(seg_path),
            Path(path),
            timetable,
            timetable_sequence,
        )
        for path in mri_paths
    ]
    pd.concat(dataframes).to_csv(output, sep=";", index=False)


def segstats_region(seg_mri, description, labels):
    region_mask = np.isin(seg_mri.data, labels)
    voxelcount = region_mask.sum()
    volscale = voxel_count_to_ml_scale(seg_mri.affine)
    return {
        "label": ",".join([str(x) for x in labels]),
        "description": description,
        "voxelcount": voxelcount,
        "volume_ml": volscale * voxelcount,
    }


def create_dataframe(
    seg_path: Path,
    mri_path: Path,
    timestamp_path: Optional[str | Path] = None,
    timestamp_sequence: Optional[str | Path] = None,
) -> pd.DataFrame:
    data_mri = sm.load_mri(mri_path, dtype=np.single)
    seg_mri = sm.load_mri(seg_path, dtype=np.int16)
    sm.assert_same_space(seg_mri, data_mri)

    seg_pattern = (
        r"(?P<subject>sub-(control|patient)*\d{2})_seg-(?P<segmentation>[^\.]+)"
    )
    mri_data_pattern = r"(?P<subject>sub-(control|patient)*\d{2})_(?P<session>ses-\d{2})_(?P<mri_data>[^\.]+)"
    lut_path = with_suffix(Path(seg_path), "_LUT.txt")
    if Path(lut_path).exists():
        lut = segtools.read_lut(lut_path)
    else:
        lut = segtools.read_lut(None)

    seg, data = seg_mri.data, data_mri.data
    seg_labels = np.unique(seg[seg != 0])
    lut_regions = lut.loc[lut.label.isin(seg_labels), ["label", "description"]].to_dict(
        "records"
    )
    regions = {
        **{d["description"]: sorted([d["label"]]) for d in lut_regions},
        **default_segmentation_groups(),
    }
    seg_info = (
        m.groupdict()
        if (m := re.match(seg_pattern, Path(seg_path).name)) is not None
        else {"segmentation": None, "subject": None}
    )
    data_info = (
        m.groupdict()
        if (m := re.match(mri_data_pattern, Path(mri_path).name)) is not None
        else {"mri_data": None, "subject": None, "session": None}
    )
    try:
        data_info["timestamp"] = find_timestamp(
            Path(str(timestamp_path)),
            str(timestamp_sequence),
            str(data_info["subject"]),
            str(data_info["session"]),
        )
    except (ValueError, RuntimeError, KeyError):
        data_info["timestamp"] = None

    info = seg_info | data_info

    records = []
    finite_mask = np.isfinite(data)
    volscale = voxel_count_to_ml_scale(seg_mri.affine)
    for description, labels in tqdm.tqdm(regions.items()):
        region_mask = np.isin(seg_mri.data, labels)
        voxelcount = region_mask.sum()
        record = {
            "label": ",".join([str(x) for x in labels]),
            "description": description,
            "voxelcount": voxelcount,
            "volume_ml": volscale * voxelcount,
        }
        if voxelcount == 0:
            records.append(record)
            continue

        data_mask = region_mask * finite_mask
        region_data = data[data_mask]
        num_nan = (~np.isfinite(region_data)).sum()
        record["num_nan_values"] = num_nan
        if num_nan == voxelcount:
            records.append(record)
            continue

        stats = {
            "sum": np.sum(region_data),
            "mean": np.mean(region_data),
            "median": np.median(region_data),
            "std": np.std(region_data),
            "min": np.min(region_data),
            **{
                f"PC{pc}": np.quantile(region_data, pc / 100)
                for pc in [1, 5, 25, 75, 90, 95, 99]
            },
            "max": np.max(region_data),
        }
        records.append({**record, **stats})

    dframe = pd.DataFrame.from_records(records)
    dframe = prepend_info(
        dframe,
        segmentation=info["segmentation"],
        mri_data=info["mri_data"],
        subject=info["subject"],
        session=info["session"],
        timestamp=info["timestamp"],
    )
    return dframe


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
            "mean": np.mean,
            "median": np.median,
            "std": np.std,
            "min": lambda x: np.min(x),
            "PC1": lambda x: np.quantile(x, 0.01),
            "PC5": lambda x: np.quantile(x, 0.05),
            "PC95": lambda x: np.quantile(x, 0.95),
            "PC99": lambda x: np.quantile(x, 0.99),
            "max": lambda x: np.max(x),
        }
        records.append(group_regions)
    return pd.DataFrame.from_records(records)


def voxel_count_to_ml_scale(affine: np.ndarray):
    return 1e-3 * np.linalg.det(affine[:3, :3])


def segstats(seg: np.ndarray, lut: pd.DataFrame, volscale: float):
    lut = lut or read_lut(None)
    labels = np.unique(seg[seg > 0])
    seg_table = pd.DataFrame.from_records(
        [
            {
                "label": label,
                "description": find_label_description(label, lut),
                "voxelcount": (seg == label).sum(),
                "volume (mL)": volscale * (seg == label).sum(),
            }
            for label in labels
        ]
    )
    total = {
        "label": set(labels),
        "description": "all-regions",
        "voxelcount": (seg != 0).sum(),
        "volume (mL)": volscale * (seg != 0).sum(),
    }
    seg_table = pd.concat(
        [
            seg_table,
            pd.DataFrame.from_records([total]),
        ],
        ignore_index=True,
    )
    return seg_table


if __name__ == "__main__":
    compute_mri_stats()
