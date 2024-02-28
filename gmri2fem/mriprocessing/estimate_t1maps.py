import argparse
from pathlib import Path

import nibabel.freesurfer.mghformat as nibmgh
import numpy as np

from gmri2fem.utils import nan_filter_gaussian

# Treshold T1 values to range [TMIN, TMAX]
TMIN, TMAX = 0.2, 4.5  # seconds


def T1_estimate_from_T1w(
    S: np.ndarray,
    S0: np.ndarray,
    T1_0: np.ndarray,
    b: float = 1.48087682,
    T_min: float = TMIN,
    T_max: float = TMAX,
) -> np.ndarray:
    # Formula is only valid for limited range.
    Smin, Smax = S0 * np.exp(b * (T1_0 - T_max)), S0 * np.exp(b * (T1_0 - T_min))
    mask = np.logical_and(
        Smin < S,
        S <= Smax,
    )
    T1 = np.nan * np.zeros_like(T1_0)
    T1[mask] = T1_0[mask] - np.log(S[mask] / S0[mask]) / b
    return T1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputT1w", type=Path, required=True)
    parser.add_argument("--referenceT1w", type=Path, required=True)
    parser.add_argument("--referenceT1map", type=Path, default=None)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sigma", type=float, default=0)
    parser.add_argument("--truncate", type=int, default=4.0)

    args = parser.parse_args()
    T1w_reference = nibmgh.load(args.referenceT1w)
    T1w_reference_affine = T1w_reference.affine
    T1w_reference_data = T1w_reference.get_fdata(dtype=np.float32)

    T1map_reference = nibmgh.load(args.referenceT1map)
    T1map_reference_data = T1map_reference.get_fdata(dtype=np.float32)
    assert np.allclose(
        T1w_reference_affine, T1map_reference.affine
    ), "Affine transformations differ, are you sure the images are registered properly?"

    T1w = nibmgh.load(args.inputT1w)
    T1w_data = T1w.get_fdata(dtype=np.float32)
    assert np.allclose(
        T1w_reference_affine, T1w.affine
    ), "Affine transformations differ, are you sure the images are registered properly?"

    T1map_data = T1_estimate_from_T1w(
        S=T1w_data, S0=T1w_reference_data, T1_0=T1map_reference_data
    )

    if args.sigma > 0:
        T1map_data = np.where(
            np.isnan(T1map_data),
            nan_filter_gaussian(T1map_data, args.sigma, args.truncate),
            T1map_data,
        )

    nibmgh.save(nibmgh.MGHImage(T1map_data, T1w_reference_affine), args.output)
