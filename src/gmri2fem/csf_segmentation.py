import os
from pathlib import Path

import numpy as np
import simple_mri as sm
import pandas as pd

from gmri2fem.masking import create_csf_mask
from gmri2fem import segment_tools as segtools
from gmri2fem.segmentation_refinement import (
    extrapolate_segmentation_to_mask,
    resample_segmentation,
)


PRECSF_EXCLUDED_SEGMENTS = [
    "Right-choroid-plexus",
    "Left-choroid-plexus",
    "WM-hypointensities",
    "CC_Posterior",
    "CC_Mid_Posterior",
    "CC_Central",
    "CC_Mid_Anterior",
    "CC_Anterior",
    "Left-Cerebral-White-Matter",
    "Right-Cerebral-White-Matter",
    "Left-Cerebellum-White-Matter",
    "Right-Cerebellum-White-Matter",
    "Left-Thalamus",
    "Right-Thalamus",
    "Left-Caudate",
    "Right-Caudate",
    "Left-Pallidum",
    "Right-Pallidum",
    "Left-Amygdala",
    "Right-Amygdala",
    "CSF",
    "Left-Accumbens-area",
    "Right-Accumbens-area",
    "Left-Putamen",
    "Right-Putamen",
]


def optional_lut(lut):
    if lut is None:
        lut = Path(os.environ["FREESURFER_HOME"]) / "FreeSurferColorLUT.txt"
    return segtools.read_lut(lut)


def translate_segment_descriptions_to_labels():
    lut = Path(os.environ["FREESURFER_HOME"]) / "FreeSurferColorLUT.txt"
    lut_table = segtools.read_lut(lut)
    return [
        segtools.find_description_label(entry, lut_table)
        for entry in PRECSF_EXCLUDED_SEGMENTS
    ]


def create_precsf_seg(aparc: np.ndarray):
    precsf_excluded_labels = translate_segment_descriptions_to_labels()
    precsf_seg = aparc.copy()
    precsf_seg[np.isin(precsf_seg, precsf_excluded_labels)] = 0
    return precsf_seg


def extend_freesurfer_lut(csf_seg: np.ndarray, fs_lut=None):
    # Create a color map adhering to FreeSurfer LUT description
    if fs_lut is None:
        freesurfer_lut = Path(os.environ["FREESURFER_HOME"]) / "FreeSurferColorLUT.txt"
    else:
        freesurfer_lut = Path(fs_lut)
    lut_table = segtools.read_lut(freesurfer_lut)

    csf_seg_labels = np.unique(csf_seg)[1:]  # Exclude 0
    csf_seg_lut = lut_table.loc[np.isin(lut_table["label"], csf_seg_labels)]
    csf_seg_lut.loc[:, "label"] += 15000
    csf_seg_lut.loc[:, "description"] = "csf-" + csf_seg_lut.loc[:, "description"]
    csf_seg_lut.loc[:, ["R", "G", "B"]] *= 0.8
    return pd.concat((lut_table, csf_seg_lut))


def create_csf_seg(aparc_mri: sm.SimpleMRI, t2w_mri: sm.SimpleMRI):
    precsf_seg = create_precsf_seg(aparc_mri.data)
    precsf_seg_mri = sm.SimpleMRI(precsf_seg, aparc_mri.affine)
    csf_mask = create_csf_mask(t2w_mri.data)

    mask_mri = sm.SimpleMRI(csf_mask, t2w_mri.affine)
    return extrapolate_segmentation_to_mask(precsf_seg_mri, mask_mri)


def create_extended_segmentation(
    aparc_mri: sm.SimpleMRI,
    t2w_mri: sm.SimpleMRI,
    wmparc_mri: sm.SimpleMRI,
):
    csf_seg_mri = create_csf_seg(aparc_mri, t2w_mri)
    csf_seg = csf_seg_mri.data

    # Both FastSurfer and FreeSurfer often resample images so that
    # they do not match the reference images space.
    wmparc_mri = resample_segmentation(wmparc_mri, csf_seg_mri)
    wmparc = wmparc_mri.data
    extended_seg = np.where(csf_seg > 0, csf_seg + 15000, wmparc)
    extended_lut = extend_freesurfer_lut(csf_seg)
    return {
        "segmentation": sm.SimpleMRI(extended_seg, t2w_mri.affine),
        "LUT": extended_lut,
    }


def create_extended_segmentation_cli(
    aparc_path: str, t2w_path: str, wmparc_path: str, output: str
):
    aparc_mri = sm.load_mri(aparc_path, dtype=np.int16)
    t2w_mri = sm.load_mri(t2w_path, dtype=np.single)
    wmparc_mri = sm.load_mri(wmparc_path, dtype=np.int16)

    extended_segmentation_d = create_extended_segmentation(
        aparc_mri, t2w_mri, wmparc_mri
    )
    extended_segmentation_mri = extended_segmentation_d["segmentation"]
    extended_lut = extended_segmentation_d["LUT"]

    sm.save_mri(extended_segmentation_mri, Path(output), dtype=np.single)
    segtools.write_lut(
        str(output).replace(".nii.gz", "_LUT.txt"),
        extended_lut,
    )

    # Also add to FS luts-directory for easy use in freeview.
    fs_base_path = (
        Path(os.environ["FREESURFER_HOME"]) / "luts/ExtendedFreeSurferColorLUT.txt"
    )
    if not fs_base_path.exists():
        segtools.write_lut(
            fs_base_path,
            extended_lut,
        )
