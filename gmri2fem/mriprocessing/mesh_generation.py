from pathlib import Path

import SVMTK as svmtk
from loguru import logger
from pantarei import mesh2xdmf, xdmf2hdf


def create_brain_mesh(stls, output, resolution, remove_ventricles=True):
    """Adjusted from
    https://github.com/kent-and/mri2fem/blob/master/mri2fem/mri2fem/chp4/fullbrain-five-domain.py
    """
    logger.info(f"Creating brain mesh from surfaces {stls}")
    surfaces = [svmtk.Surface(str(stl)) for stl in stls]

    # Merge lh rh white surface, and drop the latter.
    surfaces[2].union(surfaces[3])
    surfaces.pop(3)

    # Define identifying tags for the different regions
    tags = {"pial": 1, "white": 2, "ventricle": 3}

    smap = svmtk.SubdomainMap()
    smap.add("1000", tags["pial"])
    smap.add("0100", tags["pial"])
    smap.add("1010", tags["white"])
    smap.add("0110", tags["white"])
    smap.add("1110", tags["white"])
    smap.add("1011", tags["ventricle"])
    smap.add("0111", tags["ventricle"])
    smap.add("1111", tags["ventricle"])

    domain = svmtk.Domain(surfaces, smap)
    domain.create_mesh(resolution)
    if remove_ventricles:
        domain.remove_subdomain(tags["ventricle"])

    domain.save(str(output.with_suffix(".mesh")))
    xdmfdir = output.parent / "mesh_xdmfs"
    xdmfdir.mkdir(exist_ok=True)
    mesh2xdmf(output.with_suffix(".mesh"), xdmfdir, dim=3)
    return xdmf2hdf(xdmfdir, output)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--surfaces", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resolution", type=int, required=True)
    args = parser.parse_args()

    surfaces = ["lh.pial", "rh.pial", "lh.white", "rh.white", "ventricles"]
    stls = [Path(stl) for stl in args.surfaces]
    stls = sorted(stls, key=lambda x: surfaces.index(x.stem))
    meshfile = create_brain_mesh(
        stls=stls,
        output=Path(args.output),
        resolution=args.resolution,
        remove_ventricles=True,
    )
    logger.info(f"Generated mesh in file {meshfile}")
