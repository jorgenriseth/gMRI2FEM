from typing import Optional

import numpy as np


def orient_and_slice_image(
    data: np.ndarray, orientation: str, slice_idx: int, cutoff: Optional[float] = None
):
    if orientation == "saggital":
        plane = np.s_[slice_idx, :, :]
        im = data[plane]
    elif orientation == "transversal":
        plane = np.s_[:, slice_idx, :]
        im = np.rot90(data[plane])
    elif orientation == "coronal":
        plane = np.s_[:, :, slice_idx]
        im = data[plane].T
    else:
        raise ValueError(f"Wrong orientation, got {orientation}")
    if cutoff is not None:
        im[abs(im) <= cutoff] = np.nan
    return im
