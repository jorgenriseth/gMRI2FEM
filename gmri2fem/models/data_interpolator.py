from pathlib import Path

import dolfin as df
import numpy as np
import pantarei.interpolator as pri
from pantarei.fenicsstorage import FenicsStorage, delete_dataset


def interpolate_from_file(filepath, name, t):
    store = FenicsStorage(filepath, "r")
    tvec = store.read_timevector(name)
    bin = np.digitize(t, tvec) - 1
    C = [store.read_function(name, idx=i) for i in range(tvec.size)[bin : bin + 2]]
    interpolator = pri.vectordata_interpolator(C, tvec[bin : bin + 2])
    u = df.Function(C[0].function_space())
    u.vector()[:] = interpolator(t)
    store.close()
    return u


if __name__ == "__main__":
    # TODO: Move to test.
    """This part mainly verifies that the interpolation of the data
    works as expected. This entails
    1. Interpolating the data at a regular set of timepoint, not necessarily
        adhering to the timepoints for the data.
    2. Interpolating between the regular timepoints from above, and evvaluating
        them at the timepoints of the data. This is done directly from file."""
    some_subject = next(Path("DATA/").iterdir())
    filepath = some_subject / "FENICS/cdata_32.hdf"
    store = FenicsStorage(filepath, "a")
    tvec = store.read_timevector("cdata")
    C = [store.read_function("cdata", idx=i) for i in range(tvec.size)]
    delete_dataset(store, "/cinterp")
    interpolator = pri.vectordata_interpolator(C, tvec)
    V = C[0].function_space()
    u = df.Function(V)
    ti = 0.0
    dt = tvec[-1] / 10
    while ti < tvec[-1]:
        u.vector()[:] = interpolator(ti)  # , V)
        store.write_checkpoint(u, "cinterp", ti)
        ti += dt
    store.close()

    from gmri2fem.mriprocessing.mri2fenics import fenicsstorage2xdmf

    fenicsstorage2xdmf(store.filepath, "cinterp", "cinterp")

    store = FenicsStorage(filepath, "a")
    delete_dataset(store, "/ccheckpoints")
    for ti in tvec:
        u = interpolate_from_file(filepath, "cinterp", ti)
        store.write_checkpoint(u, "ccheckpoints", ti)
    fenicsstorage2xdmf(store.filepath, "ccheckpoints", "ccheckpoints")
    store.close()
