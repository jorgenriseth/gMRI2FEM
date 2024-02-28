import argparse
from pathlib import Path

# parser = argparse.ArgumentParser()
# parser.add_argument("--resolution", type=int, required=True)
# args = parser.parse_args()
from types import SimpleNamespace
from typing import Literal, Optional

import matplotlib.cm as mplcmaps
import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt
import nibabel
import numpy as np
from scipy.spatial import distance

from gmri2fem.mriprocessing.overlay_mri_image import orient_and_slice_image

args = SimpleNamespace(resolution=32)


def rel_voxelwise_err(D_p, D_d, f=None):
    if f is None:
        f = lambda data, target: target
    return (D_p - D_d) / f(D_p, D_d)


def rel_voxelwise_percent_difference(D_p, D_d):
    return (D_p - D_d) / ((np.abs(D_p) + np.abs(D_d)) / 2)


def median_relative_voxelwise_err(D_p, D_d):
    return (D_p - D_d) / np.nanmedian(np.abs(D_d))


def std_relative_voxelwise_err(D_p, D_d):
    return (D_p - D_d) / np.nanstd(D_d)


def max_relative_voxelwise_err(D_p, D_d):
    return (D_p - D_d) / np.nanmax(np.abs(D_d))


def global_p_norm_error(D_p, D_d, p):
    mask = ~np.isnan(D_p)
    vec_p, vec_d = D_p[mask], D_d[mask]
    return np.linalg.norm(vec_p - vec_d, ord=p) / (
        np.linalg.norm(vec_p, ord=p) * np.linalg.norm(vec_d, ord=p)
    )


def global_cosine_error(D_p, D_d):
    mask = ~np.isnan(D_p)
    vec_p, vec_d = D_p[mask], D_d[mask]
    return distance.cosine(vec_p, vec_d)


def angular_relative_difference(D_p, D_d):
    return 2 * np.arctan(D_p / D_d) - np.pi / 2


def global_score_error(D_p, D_d, f):
    mask = ~np.isnan(D_p)
    vec_p, vec_d = D_p[mask], D_d[mask]
    return (f(vec_p) - f(vec_d)) / f(vec_d)


# %% Load MRI-data
res = args.resolution
mri_raw_dir = Path("DATA/CTRL_001/CONCENTRATIONS")
mri_processed_dir = Path(f"DATA/CTRL_001/MODELING/resolution{res}/MRIs")
mri_raw_paths = sorted([x for x in mri_raw_dir.glob("*_*.mgz")])
mri_processed_paths = sorted([x for x in mri_processed_dir.glob("data*.nii.gz")])

mris_raw = [nibabel.load(x).get_fdata() for x in mri_raw_paths[1:]]
mris_processed = [nibabel.load(x).get_fdata() for x in mri_processed_paths[1:]]
mask = np.prod(
    [
        ~np.isnan(raw * processed) * ~np.isinf(raw * processed)
        for raw, processed in zip(mris_raw, mris_processed)
    ],
    axis=0,
).astype(bool)

for raw, processed in zip(mris_raw, mris_processed):
    processed[~mask] = np.nan
    raw[~mask] = np.nan

e = [proc - raw for (proc, raw) in zip(mris_processed, mris_raw)]
e_rel = [e / np.abs(raw) for e, raw in zip(e, mris_raw)]
e_median = [e / np.nanmedian(np.abs(raw)) for e, raw in zip(e, mris_raw)]
e_mean = [e / np.nanmean(np.abs(raw)) for e, raw in zip(e, mris_raw)]
e_max = [e / np.nanmax(np.abs(raw)) for e, raw in zip(e, mris_raw)]
e_rpd = [e / (0.5 * (np.abs(e + raw) + np.abs(raw))) for e, raw in zip(e, mris_raw)]

# %% Plot MRI-format comparison
raw_images = [orient_and_slice_image(mri, "transversal", 100) for mri in mris_raw]
processed_images = [
    orient_and_slice_image(mri, "transversal", 100) for mri in mris_processed
]
diff_images = [orient_and_slice_image(mri, "transversal", 100) for mri in e]

# %%
color_quantiles = 0.01
mri_color_min = np.nanquantile([mris_raw, mris_processed], color_quantiles)
mri_color_max = np.nanquantile([mris_raw, mris_processed], 1 - color_quantiles)

print(
    "MRI min, max:",
    np.nanmin([mris_raw, mris_processed]),
    np.nanmax([mris_raw, mris_processed]),
)
print(
    f"MRI-color range, quantiles {color_quantiles, 1-color_quantiles}:",
    mri_color_min,
    mri_color_max,
)

diff_quantile = 0.99
diff_color_range = np.nanquantile(np.abs(diff_images), diff_quantile)
print(
    "Error-range min, max:",
    np.nanmin(diff_images),
    np.nanmax(diff_images),
)
diff_color_range = np.nanquantile(np.abs(diff_images), diff_quantile)
print(
    f"Error-color range, abs-value quantile {diff_quantile}:",
    +diff_color_range,
)

nrows = 3
ncols = 4
cbar_width = 100
sizes = [x.shape for x in raw_images]
total_width = sum([min(x[0], x[1]) for x in sizes]) + cbar_width
total_height = nrows * max([max(x[0], x[1]) for x in sizes])
fig_width = 12
scale_factor = fig_width / total_width
fig_height = scale_factor * total_height
horizontal_offsets = [0, *np.cumsum([max(x[0], x[1]) for x in sizes]) / total_width]
vertical_offsets = 1.0 - np.cumsum([max(x[0], x[1]) for x in sizes]) / total_height
image_widths = np.diff(horizontal_offsets)
image_heights = -np.diff(vertical_offsets)
fig = plt.figure(figsize=(fig_width, fig_height))
for idx, images in enumerate(zip(raw_images, processed_images, diff_images)):
    for row_idx, im in enumerate(images):
        rect = [
            horizontal_offsets[idx],
            vertical_offsets[row_idx],
            image_widths[idx],
            image_heights[row_idx],
        ]
        ax = fig.add_axes(rect=rect)
        if row_idx == 2:
            ax.imshow(
                im, cmap="coolwarm", vmin=-diff_color_range, vmax=diff_color_range
            )
        else:
            ax.imshow(im, cmap="magma", vmin=mri_color_min, vmax=mri_color_max)
        ax.set_xticks([])
        ax.set_yticks([])
cmap = mplcmaps.magma
cax1 = fig.add_axes(
    [
        horizontal_offsets[-1],
        vertical_offsets[1],
        (1 - horizontal_offsets[-1]) / 2,
        1 - vertical_offsets[1],
    ]
)
norm = mplcolors.Normalize(vmin=mri_color_min, vmax=mri_color_max)
c = fig.colorbar(
    mplcmaps.ScalarMappable(norm=norm, cmap=cmap),
    cax=cax1,
    orientation="vertical",
    label="Concentration [mM]",
)
c.set_label(label="", size=40)
cax1.tick_params(axis="y", labelsize=16)
cmap = mplcmaps.coolwarm
cax2 = fig.add_axes(
    [
        horizontal_offsets[-1],
        vertical_offsets[2] + 0.02,
        (1 - horizontal_offsets[-1]) / 2,
        image_heights[2] - 0.04,
    ],
    xticks=[],
    yticks=[],
)
cax2.set_frame_on(False)
norm = mplcolors.Normalize(vmin=-diff_color_range, vmax=diff_color_range)
c2 = fig.colorbar(
    mplcmaps.ScalarMappable(norm=norm, cmap=cmap),
    ax=cax2,
    orientation="vertical",
    label="Relative difference",
)
c.set_label(label="", size=40)
cax2.tick_params(axis="y", labelsize=16)
plt.savefig(f"figures/mri-interpolation/diffs-resolution{res}.png")
plt.show()


# %%
errors = {
    "rel_1norm": [
        global_p_norm_error(proc, raw, 1) for proc, raw in zip(mris_processed, mris_raw)
    ],
    "rel_2norm": [
        global_p_norm_error(proc, raw, 2) for proc, raw in zip(mris_processed, mris_raw)
    ],
    "cosine_error": [
        global_cosine_error(proc, raw) for proc, raw in zip(mris_processed, mris_raw)
    ],
    "rel_mean": [
        global_score_error(proc, raw, np.nanmean)
        for proc, raw in zip(mris_processed, mris_raw)
    ],
}

plt.figure()
plt.plot(range(1, 5), errors["rel_1norm"], label="Relative 1-norm")
plt.plot(range(1, 5), errors["rel_2norm"], label="Relative 2-norm")
plt.plot(range(1, 5), errors["cosine_error"], label="Cosine Distance")
plt.legend()
plt.ylim(0, None)
plt.savefig(f"figures/mri-interpolation/errornorms-resolution{res}.png")

plt.figure()
for idx, data in enumerate(e_max):
    plt.violinplot(data[~np.isnan(data)], [idx], showmeans=True)
plt.savefig(f"figures/mri-interpolation/diff-distribution-resolution{res}")

plt.show()

# %% Plot quantities of interest

import dolfin as df
import pantarei as pr

plt.figure()
plt.plot([np.nanmean(mri_d) for mri_d in mris_raw], label="Raw")
plt.plot([np.nanmean(mri_p) for mri_p in mris_processed], label="Processed")
plt.title("Means")
plt.legend()

plt.figure()
plt.plot([np.nanmedian(mri_d) for mri_d in mris_raw], label="Raw")
plt.plot([np.nanmedian(mri_p) for mri_p in mris_processed], label="Processed")
plt.title("Medians")
plt.legend()
plt.show()


SUBDOMAIN_LABELS = {
    "gray": 1,
    "white": 2,
}


def fem_solute_quantities(
    input: Path, subdomain_labels: dict[str, int]
) -> dict[str, dict[str, int | list[int] | float | list[float]]]:
    with df.HDF5File(df.MPI.comm_world, str(input), "r") as hdf:
        domain = pr.read_domain(hdf)
        dx = df.Measure("dx", domain, subdomain_data=domain.subdomains)
        timevec = pr.read_timevector(hdf, "total_concentration")
        u = pr.read_function(hdf, "total_concentration", domain)
        region_data = {
            region: {
                "idx": region_idx,
                "volume": df.assemble(1 * dx(region_idx)),
                "mass": np.nan * np.ones_like(timevec),
                "mean": np.nan * np.ones_like(timevec),
            }
            for region, region_idx in subdomain_labels.items()
        }
        for idx, ti in enumerate(timevec):
            pr.read_checkpoint(hdf, u, "total_concentration", idx)
            for region, data in region_data.items():
                data["mass"][idx] = df.assemble(u * dx(data["idx"]))
                data["mean"][idx] = data["mass"][idx] / data["volume"]

    region_data["whole"] = {
        "idx": list(data["idx"] for data in region_data.values()),
        "volume": sum((data["volume"] for data in region_data.values())),
        "mass": sum((data["mass"] for data in region_data.values())),
    }
    region_data["whole"]["mean"] = (
        region_data["whole"]["mass"] / region_data["whole"]["volume"]
    )
    return region_data


hdffile = Path(f"DATA/CTRL_001/MODELING/resolution{res}/data.hdf")
region_fem_data = fem_solute_quantities(hdffile, SUBDOMAIN_LABELS)

# %% Regionwise
from gmri2fem.analysis.seg_groups import default_segmentation_groups

asegfile = "DATA/CTRL_001/freesurfer/mri/aseg.mgz"
aseg = nibabel.load(asegfile).get_fdata().astype(int)
groups = default_segmentation_groups()
group_masks = {
    "white": np.isin(aseg, groups["white-matter"]),
    "gray": np.isin(aseg, groups["gray-matter"]),
    "whole": mask,
}

fig, axes = plt.subplots(3, 1, figsize=(12, 6))
fig.suptitle("Means")
axes[0].plot(
    [np.nanmean(mri_d[group_masks["white"]]) for mri_d in mris_raw], label="Raw"
)
axes[0].plot(
    [np.nanmean(mri_p[group_masks["white"]]) for mri_p in mris_processed],
    label="Processed",
)
axes[0].plot(region_fem_data["white"]["mean"][1:], label="FEM-format")
axes[0].set_ylabel("White")
axes[0].legend()

axes[1].plot(
    [np.nanmean(mri_d[group_masks["gray"]]) for mri_d in mris_raw], label="Raw"
)
axes[1].plot(
    [np.nanmean(mri_p[group_masks["gray"]]) for mri_p in mris_processed],
    label="Processed",
)
axes[1].plot(region_fem_data["gray"]["mean"][1:], label="FEM-format")
axes[1].set_ylabel("Gray")
axes[1].legend()

axes[2].plot([np.nanmean(mri_d) for mri_d in mris_raw], label="Raw")
axes[2].plot([np.nanmean(mri_p) for mri_p in mris_processed], label="Processed")
axes[2].plot(region_fem_data["whole"]["mean"][1:], label="FEM-format")
axes[2].set_ylabel("Whole")
axes[2].legend()
plt.savefig(
    "figures/mri-interpolation/mean-concentration-region-compare.pdf",
    bbox_inches="tight",
)

fig, axes = plt.subplots(3, 1, figsize=(12, 6))
fig.suptitle("Mean Errors")
axes[0].plot(
    [np.nanmean(mri_d[group_masks["white"]]) for mri_d in mris_raw], label="Raw"
)
axes[0].plot(
    [np.nanmean(mri_p[group_masks["white"]]) for mri_p in mris_processed],
    label="Processed",
)
axes[0].plot(region_fem_data["white"]["mean"][1:], label="FEM-format")
axes[0].set_ylabel("White")
axes[0].legend()

axes[1].plot(
    [np.nanmean(mri_d[group_masks["gray"]]) for mri_d in mris_raw], label="Raw"
)
axes[1].plot(
    [np.nanmean(mri_p[group_masks["gray"]]) for mri_p in mris_processed],
    label="Processed",
)
axes[1].plot(region_fem_data["gray"]["mean"][1:], label="FEM-format")
axes[1].set_ylabel("Gray")
axes[1].legend()

axes[2].plot([np.nanmean(mri_d) for mri_d in mris_raw], label="Raw")
axes[2].plot([np.nanmean(mri_p) for mri_p in mris_processed], label="Processed")
axes[2].plot(region_fem_data["whole"]["mean"][1:], label="FEM-format")
axes[2].set_ylabel("Whole")
axes[2].legend()
plt.savefig(
    "figures/mri-interpolation/mean-concentration-region-compare.pdf",
    bbox_inches="tight",
)

E_white = [
    (np.nanmean(mri_p[group_masks["white"]]) - np.nanmean(mri_d[group_masks["white"]]))
    / np.nanmean(mri_d[group_masks["white"]])
    for mri_p, mri_d in zip(mris_processed, mris_raw)
]


fig, axes = plt.subplots(3, 1, figsize=(12, 6))
fig.suptitle("Mean Errors")
axes[0].plot(E_white)
axes[0].set_ylabel("White")
axes[0].legend()
plt.show()

axes[1].plot(
    [np.nanmean(mri_d[group_masks["gray"]]) for mri_d in mris_raw], label="Raw"
)
axes[1].plot(
    [np.nanmean(mri_p[group_masks["gray"]]) for mri_p in mris_processed],
    label="Processed",
)
axes[1].plot(region_fem_data["gray"]["mean"][1:], label="FEM-format")
axes[1].set_ylabel("Gray")
axes[1].legend()

axes[2].plot([np.nanmean(mri_d) for mri_d in mris_raw], label="Raw")
axes[2].plot([np.nanmean(mri_p) for mri_p in mris_processed], label="Processed")
axes[2].plot(region_fem_data["whole"]["mean"][1:], label="FEM-format")
axes[2].set_ylabel("Whole")
axes[2].legend()
plt.savefig(
    "figures/mri-interpolation/mean-concentration-region-compare.pdf",
    bbox_inches="tight",
)


fig, axes = plt.subplots(3, 1, figsize=(12, 6))
fig.suptitle("Medians")
axes[0].plot(
    [np.nanmedian(mri_d[group_masks["white"]]) for mri_d in mris_raw], label="Raw"
)
axes[0].plot(
    [np.nanmedian(mri_p[group_masks["white"]]) for mri_p in mris_processed],
    label="Processed",
)
axes[0].set_ylabel("White")
axes[0].legend()

axes[1].plot(
    [np.nanmedian(mri_d[group_masks["gray"]]) for mri_d in mris_raw], label="Raw"
)
axes[1].plot(
    [np.nanmedian(mri_p[group_masks["gray"]]) for mri_p in mris_processed],
    label="Processed",
)
axes[1].set_ylabel("Gray")
axes[1].legend()

axes[2].plot([np.nanmedian(mri_d) for mri_d in mris_raw], label="Raw")
axes[2].plot([np.nanmedian(mri_p) for mri_p in mris_processed], label="Processed")
axes[2].set_ylabel("Whole")
axes[2].legend()
plt.savefig(
    "figures/mri-interpolation/median-concentration-region-compare.pdf",
    bbox_inches="tight",
)
plt.show()
