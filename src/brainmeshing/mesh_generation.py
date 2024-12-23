from pathlib import Path

import click
import SVMTK as svmtk
from pantarei.meshprocessing import mesh2xdmf, xdmf2hdf

from brainmeshing.utils import subdomain_mapper


@click.command("meshgen")
@click.option("--surfacedir", type=Path, required=True)
@click.option("--output", type=Path, required=True)
@click.option("--resolution", type=int, required=True)
@click.option("--keep-ventricles", is_flag=True)
def meshgen(surfacedir: Path, output: Path, resolution: float, keep_ventricles: bool):
    surface_files = {
        surf: surfacedir / f"{surf}.stl"
        for surf in [
            "rh_pial_novent",
            "lh_pial_novent",
            "subcortical_gm",
            "white",
            "ventricles",
        ]
    }
    for surf in surface_files.values():
        assert surf.exists(), f"Missing surface file, {surf}"
    try:
        svmtk_surfaces = {
            surf: svmtk.Surface(str(path)) for surf, path in surface_files.items()
        }
    except Exception as e:
        print(surface_files)
        raise e
    tags = {"gray": 1, "white": 2, "subcort-gray": 3, "ventricles": 4}
    surfaces = [
        svmtk_surfaces[surf]
        for surf in [
            "rh_pial_novent",
            "lh_pial_novent",
            "subcortical_gm",
            "white",
            "ventricles",
        ]
    ]

    smap = svmtk.SubdomainMap(num_surfaces=len(surfaces))
    subdomain_mapper(smap, "....1", tags["ventricles"])
    subdomain_mapper(smap, "..1.0", tags["subcort-gray"])
    subdomain_mapper(smap, "..010", tags["white"])
    subdomain_mapper(smap, "1.000", tags["gray"])
    subdomain_mapper(smap, "01000", tags["gray"])

    domain = svmtk.Domain(surfaces, smap)
    domain.create_mesh(resolution)
    if not keep_ventricles:
        domain.remove_subdomain(tags["ventricles"])
    domain.save(str(output.with_suffix(".mesh")))

    xdmfdir = output.parent / "mesh_xdmfs"
    xdmfdir.mkdir(exist_ok=True, parents=True)
    mesh2xdmf(output.with_suffix(".mesh"), xdmfdir, dim=3)
    hdf_out = xdmf2hdf(xdmfdir, output)
    print(hdf_out)


if __name__ == "__main__":
    meshgen()
