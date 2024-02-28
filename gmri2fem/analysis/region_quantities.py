import json
from pathlib import Path

import nibabel
import numpy as np
import pandas as pd

from gmri2fem.analysis.seg_groups import default_segmentation_groups
from gmri2fem.filters import is_T1_mgz

# Intervals used to group the varying timestamps of different subjects
# together. Intervals are given in hours, and data gets the label 't{i}'
# i.e. t=4.5h -> 't1', t=23.2h -> 't2', and so on.
SESSION_INTERVALS = [(0, 0), (2, 6), (20, 28), (45, 52), (68, 72)]


def timestamp_session_idx(t: float) -> int:
    """Get the index of the session corresponding to the given timestamp."""
    for idx, interval in enumerate(SESSION_INTERVALS):
        if interval[0] <= t / 3600 <= interval[1]:
            return idx
    raise ValueError(f"Timestamp {t}s={t/3600:.2f}h is not in any session intervals")


def build_lut_dataframe(
    asegfile: Path, concentrationdir: Path, lut_path: Path, timestampfile: Path
):
    with open(lut_path, "r") as f:
        fs_lut = json.load(f)

    # Extract all freesurfer-lookuptable values present in aseg
    aseg = nibabel.load(asegfile).get_fdata().astype(int)
    aseg_labels = np.unique(aseg[aseg != 0])
    regions = [fs_lut[str(label)] for label in aseg_labels]
    counts = [len(aseg[aseg == label]) for label in aseg_labels]
    lut_values = pd.Series(aseg_labels)
    lut_regions = pd.Series(regions)
    lut_size = pd.Series(counts)
    timestamp = np.loadtxt(timestampfile)

    data_series = {}
    concentrations = sorted(filter(is_T1_mgz, concentrationdir.iterdir()))
    for i, cpath in enumerate(concentrations):
        concentration_image = nibabel.load(cpath).get_fdata()

        contents = np.zeros_like(aseg_labels, dtype=float)
        means = np.nan * np.zeros_like(aseg_labels, dtype=float)
        medians = np.nan * np.zeros_like(aseg_labels, dtype=float)
        stds = np.nan * np.zeros_like(aseg_labels, dtype=float)
        for idx, label in enumerate(aseg_labels):
            label_concentrations = concentration_image[aseg == label]
            label_concentrations = label_concentrations[~np.isnan(label_concentrations)]
            contents[idx] = np.sum(label_concentrations)
            if label_concentrations.size == 0:
                continue
            means[idx] = np.mean(label_concentrations)
            medians[idx] = np.median(label_concentrations)
            stds[idx] = np.std(label_concentrations)
        session_idx = timestamp_session_idx(timestamp[i])
        data_series[f"t{session_idx}"] = timestamp[i]
        data_series[f"content-t{session_idx}"] = pd.Series(contents)
        data_series[f"mean-t{session_idx}"] = pd.Series(means)
        data_series[f"median-t{session_idx}"] = pd.Series(medians)
        data_series[f"std-t{session_idx}"] = pd.Series(stds)

    return pd.DataFrame(
        {
            "FS_LUT-values": lut_values,
            "FS_LUT-regions": lut_regions,
            f"FS_LUT-voxelcount": lut_size,
            **data_series,
        }
    )


def build_groups_dataframe(asegfile: Path, concentrationdir: Path, timestampfile: Path):
    aseg = nibabel.load(asegfile).get_fdata().astype(int)

    timestamp = np.loadtxt(timestampfile)
    groups = default_segmentation_groups()
    lut_values = pd.Series(groups.values())
    lut_regions = pd.Series(groups.keys())
    lut_size = pd.Series(
        [len(aseg[np.isin(aseg, values)]) for values in groups.values()]
    )
    N = len(groups)

    data_series = {}
    concentrations = sorted(filter(is_T1_mgz, concentrationdir.iterdir()))
    for i, cpath in enumerate(concentrations):
        concentration_image = nibabel.load(cpath).get_fdata()
        contents = np.zeros(N)
        means = np.nan * np.zeros(N)
        medians = np.nan * np.zeros(N)
        stds = np.nan * np.zeros(N)
        for idx, (group, values) in enumerate(groups.items()):
            group_mask = np.isin(aseg, values)
            group_concentrations = concentration_image[group_mask]
            group_concentrations = group_concentrations[~np.isnan(group_concentrations)]

            contents[idx] = np.sum(group_concentrations)
            if group_concentrations.size == 0:
                continue
            means[idx] = np.mean(group_concentrations)
            medians[idx] = np.median(group_concentrations)
            stds[idx] = np.std(group_concentrations)

        session_idx = timestamp_session_idx(timestamp[i])
        data_series[f"t{session_idx}"] = timestamp[i]
        data_series[f"content-t{session_idx}"] = pd.Series(contents)
        data_series[f"mean-t{session_idx}"] = pd.Series(means)
        data_series[f"median-t{session_idx}"] = pd.Series(medians)
        data_series[f"std-t{session_idx}"] = pd.Series(stds)

    return pd.DataFrame(
        {
            "FS_LUT-values": lut_values,
            "FS_LUT-regions": lut_regions,
            f"FS_LUT-voxelcount": lut_size,
            **data_series,
        }
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--asegfile", type=str, required=True)
    parser.add_argument("--timestamps", type=str, required=True)
    parser.add_argument("--concentrationdir", type=str, required=True)
    parser.add_argument("--lutfile", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()

    asegfile = Path(args.asegfile)
    concentrationdir = Path(args.concentrationdir)
    lut_path = Path(args.lutfile)
    timestampfile = Path(args.timestamps)

    dframe1 = build_lut_dataframe(asegfile, concentrationdir, lut_path, timestampfile)
    dframe2 = build_groups_dataframe(asegfile, concentrationdir, timestampfile)
    dframe = pd.concat([dframe1, dframe2], ignore_index=True)
    dframe = dframe.reindex(sorted(dframe.columns), axis=1)
    dframe.to_csv(Path(args.output), index=False)
