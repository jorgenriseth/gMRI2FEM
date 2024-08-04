import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import nibabel
import numpy as np
import pydicom
import scipy


def estimate_T1_R1(filename: Path, TRESHOLDCSF: int = 220):
    # Read dicom and metadata
    volume, info = read_enhanced_dicom_volume(filename)

    # Extract some initial information
    n_slices = info[0x2001, 0x1018].value
    _, TRSE = info[0x2005, 0x1030].value
    TI = info.InversionTime
    TE = info.PerFrameFunctionalGroupsSequence[0][0x0018, 0x9112][0][0x0018, 0x0080].value  # fmt: skip

    # Create look-up table for signal-fraction to T1-estimate.
    fractionCurve, t1_grid = calculateFractionCurve(TRSE, TI, TE)

    # # Extract spin echo modulus images and inversion recovery phase corrected real images
    # to allow for phase insensitive T1 estimates
    imgVolumeSE = volume[3 * n_slices : 4 * n_slices, :, :]
    imgVolumeIR = volume[2 * n_slices : 3 * n_slices, :, :]

    # Calculating fraction
    fractionImgVolume = imgVolumeSE - imgVolumeIR / (imgVolumeSE + imgVolumeIR)

    # Estimate T1-times from interpolating faction image from look-up table
    interpolator = scipy.interpolate.interp1d(
        t1_grid, fractionCurve, kind="nearest", fill_value="extrapolate"
    )
    estimated_T1_milliseconds = interpolator(fractionImgVolume).astype(np.float32)

    # Calulate the relaxation rate (1/ms)
    estimated_R1_milliseconds = 1 / estimated_T1_milliseconds

    # Mask noisy picels based on treshold in SE image
    TRESHOLDCSF = 220
    estimated_T1_milliseconds[imgVolumeSE < TRESHOLDCSF] = 0
    estimated_R1_milliseconds[imgVolumeSE < TRESHOLDCSF] = 0

    return estimated_T1_milliseconds, estimated_R1_milliseconds


def read_enhanced_dicom_volume(filename: Path) -> tuple[np.ndarray, pydicom.Dataset]:
    info = pydicom.dcmread(filename)
    num_frames = int(info.NumberOfFrames)
    pixel_data = info.pixel_array.astype(np.single)
    volume = np.zeros_like(pixel_data)

    # For each frame, scale pixel-data with slope and intercept
    per_frame_groups = info.PerFrameFunctionalGroupsSequence
    slopes = set()
    for idx in range(num_frames):
        pfg_i = per_frame_groups[idx]
        real_world_mapping_sequence = pfg_i[0x0040, 0x9096][0]
        slope = np.float32(real_world_mapping_sequence[0x0040, 0x9225].value)
        intercept = np.float32(real_world_mapping_sequence[0x0040, 0x9224].value)
        volume[idx, :, :] = pixel_data[idx, :, :] * slope + intercept
        if slope not in slopes:
            slopes.add(slope)
            print(idx, slope, intercept)
    # volume = pixel_data.copy()

    return volume, info


def calculateFractionCurve(TRSE, TI, TE):
    """Calulate a lookup table relating a measured signal fraction
    to a T1-estimate

    Tables sequence timing parameters of the mixed sequence as input
    and caclulates the fraction
        (sSE-sIR)/(sSE+SIR)
    for T1-values between LOWTILIMIT and HIGHTLIMIT with 1ms stepping.

    The mixed sequence is a combination of two MRI acquisitions which
    run interleaced so that two separate contrasts are obtained with
    different T1-weighting The measured signales may be estimated by
    the signal equations:

    SSE = M0 * (1 - exp(-TR / T1))

    for the spin echo signal.
    """
    # Create look-up table for signal-fraction to T1-estimate.
    LOWT1LIMIT = 100
    HIGHT1LIMIT = 3500
    t1_grid = np.arange(LOWT1LIMIT, HIGHT1LIMIT + 1)
    M0 = 1
    nominalTR = TRSE
    readoutduration = 2 * TE
    TR = nominalTR - readoutduration
    sSE = M0 * (1 - np.exp(-TR / t1_grid))
    sIR = M0 - (M0 + sSE) * np.exp(-TI / t1_grid)
    fractionCurve = (sSE - sIR) / (sSE + sIR)
    return fractionCurve, t1_grid


def resample_image_volume(a: np.ndarray, dim_x: int, dim_y: int, dim_z: int):
    raise NotImplementedError("Resampling not yet implemented.")


def get_image_affine(dicompath: Path):
    """Get the affine header for nii-output.
    TODO: Parse directly from dicom-header instead.
    """
    tmpdir = Path(tempfile.mkdtemp())
    cmd = f"dcm2niix -f mixed -o '{tmpdir}' '{dicompath}'"
    subprocess.run(cmd, shell=True)
    nii = nibabel.nifti1.load(tmpdir / "mixed_e1.nii")
    shutil.rmtree(tmpdir)
    return nii.affine


def process_mixed_dicom(
    filepath: Path,
    t1output: Path,
    r1output: Path,
    TRESHOLDCSF: int = 220,
    resample: Optional[tuple[int, int, int]] = None,
):
    estimated_T1_milliseconds, estimated_R1_milliseconds = estimate_T1_R1(
        filename=filepath, TRESHOLDCSF=TRESHOLDCSF
    )
    affine = get_image_affine(filepath)

    # flip data axes to match niftii coordinates, and save as niftii
    # could potentially be fixed by reading affine from dicom, rather
    # than current approach with using dcm2niix
    t1data = estimated_T1_milliseconds
    t1data = np.flip(np.rot90(t1data, axes=(2, 1)), axis=1)
    t1 = nibabel.nifti1.Nifti1Image(t1data, affine)
    nibabel.nifti1.save(t1, t1output)

    r1data = estimated_R1_milliseconds
    r1data = np.flip(np.rot90(r1data, axes=(2, 1)), axis=1)
    r1 = nibabel.nifti1.Nifti1Image(r1data, affine)
    nibabel.nifti1.save(r1, r1output)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--dicom_session_path", type=Path, required=True)
    parser.add_argument("--output_t1", type=Path, required=True)
    parser.add_argument("--output_r1", type=Path, required=True)
    args = parser.parse_args()

    from process_dicom import (
        find_sequence_dicomdir,
        subject_sessions,
        generate_metadata,
    )

    subjects = generate_metadata("participants-private.json")
    sessions = subject_sessions()

    mixed_dir = find_sequence_dicomdir(
        args.dicom_session_path, r"(DICOM_\d+_\d+[_ ])(.*(Mixed).*)"
    )
    process_mixed_dicom(mixed_dir / "DICOM/IM_0002", args.output_t1, args.output_r1)

    # Read dicom and metadata

    exit()


if False:
    # %%
    dicom_session_path = "sub-01/DICOM-full/2023_02_13/Zi_121927/"
    # dicom_session_path = "sub-01/DICOM-full/2023_02_13/Zi_121304"
    dicom_session_path = Path(dicom_session_path)
    mixed_dir = find_sequence_dicomdir(
        Path(dicom_session_path), r"(DICOM_\d+_\d+[_ ])(.*(Mixed).*)"
    )

    filename = mixed_dir / "DICOM/IM_0002"
    # volume, info = read_enhanced_dicom_volume(filename)

    info = pydicom.dcmread(filename)
    num_frames = int(info.NumberOfFrames)
    pixel_data = info.pixel_array.astype(np.single)
    volume = pixel_data.copy()

    # # # For each frame, scale pixel-data with slope and intercept
    # per_frame_groups = info.PerFrameFunctionalGroupsSequence
    # slopes = set()
    # f in range() r idx in range(num_frames):
    #     pfg_i = per_frame_groups[idx]
    #     real_world_mapping_sequence = pfg_i[0x0040, 0x9096][0]
    #     slope = np.float32(real_world_mapping_sequence[0x0040, 0x9225].value)
    #     intercept = np.float32(real_world_mapping_sequence[0x0040, 0x9224].value)
    #     volume[idx, :, :] = pixel_data[idx, :, :] * slope + intercept
    #     if slope not in slopes:
    #         slopes.add(slope)
    #         print(idx, slope, intercept)
    #     # volume = pixel_data.copy()
    #
    # return volume, info
    # # Extract spin echo modulus images and inversion recovery phase corrected real images
    # to allow for phase insensitive T1 estimates

    n_slices = info[0x2001, 0x1018].value
    fig, axes = plt.subplots(2, 3)
    for i in range(6):
        vol = volume[i * n_slices : (i + 1) * n_slices]
        axes[i // 3, i % 3].imshow(vol[slice_])
    plt.show()

    imgVolumeSE = volume[3 * n_slices : 4 * n_slices, :, :]
    imgVolumeIR = volume[2 * n_slices : 3 * n_slices, :, :]

    print(imgVolumeIR.mean(), imgVolumeSE.mean())

    # %%
    # Extract some initial information
    TRIR, TRSE = info[0x2005, 0x1030].value
    TI = info.InversionTime
    TE = info.PerFrameFunctionalGroupsSequence[0][0x0018, 0x9112][0][0x0018, 0x0080].value  # fmt: skip

    # Create look-up table for signal-fraction to T1-estimate.
    fractionCurve, t1_grid = calculateFractionCurve(TRSE, TI, TE)

    # # Extract spin echo modulus images and inversion recovery phase corrected real images
    # to allow for phase insensitive T1 estimates
    imgVolumeSE = volume[3 * n_slices : 4 * n_slices, :, :]
    imgVolumeIR = volume[2 * n_slices : 3 * n_slices, :, :]

    # Calculating fraction
    fractionImgVolume = imgVolumeSE - imgVolumeIR / (imgVolumeSE + imgVolumeIR)

    import matplotlib.pyplot as plt
    import nibabel.nifti1 as nibnifti
    import numpy as np

    # from mixed_to_t1map import get_image_affine, read_enhanced_dicom_volume

    # %%
    slice_ = np.s_[:, 150, :]
    fig, axes = plt.subplots(1, 2)
    axes[0].imshow(imgVolumeSE[slice_])
    axes[1].imshow(imgVolumeIR[slice_])

    descs = ["e1", "e1_reala", "e1_real", "e1a", "e2"]
    fig, axes = plt.subplots(1, 6)
    idx = 0
    slice_idx = 270
    for desc in descs:
        slice_ = np.s_[:, :, slice_idx]
        dcm2nii_file = f"sub-01/ses-02/anat/sub-01_ses-02_Mixed_{desc}.nii"
        dcm2nii_vol = nibnifti.load(dcm2nii_file)
        D = dcm2nii_vol.get_fdata()
        print(desc, D.shape)
        try:
            axes[idx].imshow(D[slice_])
            idx += 1
        except:
            slice_ = np.s_[:, :, slice_idx]
            D1, D2 = D.transpose(3, 2, 0, 1)
            print(D1.shape, D2.shape)
            axes[idx].imshow(D1[slice_])
            axes[idx + 1].imshow(D2[slice_])
            idx += 2
    plt.show()
