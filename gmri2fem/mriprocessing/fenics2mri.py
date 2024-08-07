import argparse
import itertools
from pathlib import Path
from typing import Optional

import dolfin as df
import nibabel
import numpy as np
from pantarei import interpolate_from_file
from tqdm import tqdm

from gmri2fem.utils import apply_affine


def measure_function(
    f: df.Function | df.Expression | df.UserExpression,
    img: nibabel.freesurfer.mghformat.MGHImage,
    bbox: Optional[tuple[np.ndarray, np.ndarray]] = None,
) -> nibabel.freesurfer.mghformat.MGHImage | nibabel.nifti1.Nifti1Image:
    """Measure function from 3D-domain in a regular 3D-grid through interpolation
    of grid-centerpoints"""
    aff = img.header.get_vox2ras()
    shape = img.dataobj.shape
    inds = np.fromiter(
        itertools.product(*(np.arange(ni) for ni in shape)), dtype=np.dtype((int, 3))
    )
    Z = apply_affine(aff, inds)

    if bbox is not None:
        relevant_coords = np.logical_and(
            np.all(bbox[0] <= Z, axis=-1), np.all(Z <= bbox[1], axis=-1)
        )
        Z = Z[relevant_coords]
        inds = inds[relevant_coords]
    D = np.nan * np.zeros_like(img.dataobj)
    for (i, j, k), z in zip(tqdm(inds), Z):
        try:
            D[i, j, k] = f(z)
        except RuntimeError as e:
            if "set_allow_extrapolation" in str(e):
                D[i, j, k] = np.nan
    return nibabel.nifti1.Nifti1Image(D.astype(np.single), aff)
    # return nibabel.freesurfer.mghformat.MGHImage(D.astype(np.single), aff)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--simulationfile", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--referenceimage", type=Path, required=True)
    parser.add_argument("--timestamps", type=Path, required=True)
    parser.add_argument("--timeidx", type=int, required=True)
    parser.add_argument("--functionname", type=str, default="concentration")
    parser.add_argument("--extrapolation_value", type=float, default=float("nan"))
    args = parser.parse_args()

    # Map to reference volume
    reference_volume = nibabel.freesurfer.mghformat.load(args.referenceimage)
    timevec = np.loadtxt(args.timestamps)
    ti = timevec[args.timeidx]
    ci = interpolate_from_file(args.simulationfile, args.functionname, ti)

    # Create domain bounding-box
    domain = ci.function_space().mesh()
    coords = domain.coordinates()
    pmin = coords.min(axis=0)
    pmax = coords.max(axis=0)

    volume = measure_function(ci, reference_volume, (pmin, pmax))
    nibabel.save(volume, args.output)
