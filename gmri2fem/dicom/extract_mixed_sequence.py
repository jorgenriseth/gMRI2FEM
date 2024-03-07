import nibabel
import pydicom
import json
import numpy as np
import matplotlib.pyplot as plt
import scipy
from pathlib import Path

VOLUME_LABELS = [
    "IR-modulus",
    "IR-real",
    "IR-phase_corrected",
    "SE-modulus",
    "SE-real",
    "T1-scanner",
]


def subject_session_dirs(subject_dicomdir: Path) -> list[Path]:
    return [
        session_dir
        for date_dir in sorted(subject_dicomdir.glob("202[34]*"))
        for session_dir in sorted(date_dir.iterdir())
    ]


def frame_group_plane_position(pfg: pydicom.Dataset) -> np.ndarray:
    plane_position_ds = pfg.PlanePositionSequence[0]
    plane_position_de = plane_position_ds[0x0020, 0x0032]
    plane_position_val = plane_position_de.value
    plane_position = np.array([float(x) for x in plane_position_val])
    return plane_position


def frame_group_plane_orientation(pfg: pydicom.Dataset) -> np.ndarray:
    plane_orientation_ds = pfg.PlaneOrientationSequence[0]
    plane_orientation_de = plane_orientation_ds[0x0020, 0x0037]
    plane_orientation_val = plane_orientation_de.value
    plane_orientation = np.array([float(x) for x in plane_orientation_val])
    F = np.zeros((3, 2))
    F[:, 0] = plane_orientation[3:]
    F[:, 1] = plane_orientation[:3]
    return F


def build_ras_affine(
    per_frame_groups: pydicom.Sequence, slice_idx_start: int, slice_idx_end: int
):
    pfg = per_frame_groups[slice_idx_start]
    dr, dc = (float(x) for x in pfg.PixelMeasuresSequence[0].PixelSpacing)
    F = frame_group_plane_orientation(pfg)
    T_1 = frame_group_plane_position(pfg)
    N = slice_idx_end - slice_idx_start
    T_N = frame_group_plane_position(per_frame_groups[slice_idx_end - 1])

    A = np.zeros((4, 4))
    A[:3, 0] = F[:, 0] * dr
    A[:3, 1] = F[:, 1] * dc
    A[:3, 2] = (T_N - T_1) / (N - 1)
    A[:3, 3] = T_1
    A[3, 3] = 1.0

    return A


def read_volume_data(
    pixel_data: np.ndarray,
    per_frame_groups: pydicom.Sequence,
    slice_idx_start: int,
    slice_idx_end: int,
) -> np.ndarray:
    N = slice_idx_end - slice_idx_start
    nrows, ncols = pixel_data.shape[1], pixel_data.shape[2]
    volume = np.nan * np.zeros((nrows, ncols, N), dtype=np.single)
    pfg = per_frame_groups[slice_idx_start]
    rwv_mapping_sequence = pfg.RealWorldValueMappingSequence[0]
    slope = rwv_mapping_sequence[0x0040, 0x9225].value
    intercept = rwv_mapping_sequence[0x0040, 0x9224].value
    for idx in range(N):
        slice_idx = idx + slice_idx_start
        volume[:, :, idx] = intercept + slope * pixel_data[slice_idx]
    return volume


def extract_volume(
    info: pydicom.Dataset, volume_idx: int
) -> nibabel.nifti1.Nifti1Image:
    num_frames_total = int(info.NumberOfFrames)
    N = num_frames_total // len(VOLUME_LABELS)
    per_frame_groups = info.PerFrameFunctionalGroupsSequence
    slice_idx_start, slice_idx_end = volume_idx * N, (volume_idx + 1) * N
    volume = read_volume_data(
        info.pixel_array.astype(np.single),
        per_frame_groups,
        slice_idx_start,
        slice_idx_end,
    )
    affine = build_ras_affine(per_frame_groups, slice_idx_start, slice_idx_end)
    return nibabel.nifti1.Nifti1Image(volume, affine)


def dicom2nii_mixed(sequencedir: Path, outputdir: Path, label: str):
    mixed_paths = sorted(sequencedir.glob("DICOM/IM_*"))
    assert len(mixed_paths) == 1
    info = pydicom.dcmread(mixed_paths[0])
    for volume_idx, volume_label in enumerate(VOLUME_LABELS):
        nii = extract_volume(info, volume_idx)
        nibabel.nifti1.save(nii, outputdir / f"{label}_{volume_label}.nii.gz")


def T1_lookup_table(M0, TRse, TI, TE, T1_low, T1_hi):
    T1_grid = np.arange(T1_low, T1_hi + 1)
    TR = TRse - 2 * TE
    Sse = M0 * (1 - np.exp(-TR / T1_grid))
    Sir = M0 - (M0 + Sse) * np.exp(-TI / T1_grid)
    fractionCurve = (Sse - Sir) / (Sse + Sir)
    return fractionCurve, T1_grid


def estimate_T1_mixed(dicom_path: Path, T1_lo: float, T1_hi: float):
    info = pydicom.dcmread(dicom_path)
    per_frame_groups = info.PerFrameFunctionalGroupsSequence
    pfg = per_frame_groups[0]
    M0 = info.MagneticFieldStrength
    TR_ir, TR_se = info[0x2005, 0x1030].value
    TI = info.InversionTime
    TE = pfg[0x0018, 0x9114][0].EffectiveEchoTime
    F, T1_grid = T1_lookup_table(M0, TR_se, TI, TE, 200, 5000)
    IR_nii = extract_volume(info, VOLUME_LABELS.index("IR-phase_corrected"))
    IR = IR_nii.get_fdata(dtype=np.single)
    SE = extract_volume(info, VOLUME_LABELS.index("SE-modulus")).get_fdata(
        dtype=np.single
    )
    F_data = (SE - IR) / (SE + IR)
    interpolator = scipy.interpolate.interp1d(
        F, T1_grid, kind="nearest", bounds_error=False, fill_value=np.nan
    )
    T1_volume = interpolator(F_data).astype(np.float32)
    return nibabel.nifti1.Nifti1Image(T1_volume, IR_nii.affine)
