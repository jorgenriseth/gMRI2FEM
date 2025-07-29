import functools
import time as pytime
from pathlib import Path

import click
import numpy as np
import pyvista as pv
import simple_mri as sm
import skimage
import SVMTK as svmtk
from loguru import logger
from panta_rhei.meshprocessing import mesh2xdmf, xdmf2hdf

from brainmeshing.utils import (
    binary_image_surface_extraction,
    embed_overlapping,
    fs_surf_to_stl,
    grow_white_connective_tissue,
    pyvista2svmtk,
    repair_disconnected_triangulation,
    repair_svmtk,
    repair_triangulation,
    subdomain_mapper,
    surface_union,
    svmtk2pyvista,
    taubin_svmtk,
)
from brainmeshing.ventricles import extract_ventricle_surface
from gmri2fem.segmentation_groups import default_segmentation_groups
from gmri2fem.utils import segmentation_smoothing


@click.command("process-surfaces")
@click.option("--fs_dir", type=Path, required=True)
@click.option("--surface_dir", type=Path, required=True)
def process_surfaces(**kwargs):
    process_cerebral_surfaces(**kwargs)


@click.command("meshgen")
@click.option("--surface_dir", type=Path, required=True)
@click.option("--output", type=Path, required=True)
@click.option("--resolution", type=int, required=True)
@click.option("--keep-ventricles", is_flag=True)
def meshgen(**kwargs):
    generate_mesh(**kwargs)


@click.command("extract-ventricles")
@click.option("--fs_dir", type=Path, required=True)
@click.option("--surface_dir", type=Path, required=True)
def extract_ventricles(
    fs_dir: Path,
    surface_dir: Path,
):
    Path(surface_dir).mkdir(exist_ok=True)
    seg_mri = sm.load_mri(Path(fs_dir) / "mri/aseg.mgz", dtype=np.int16)
    ventricles = extract_ventricle_surface(
        seg_mri,
        initial_smoothing=0,
        min_radius=3,
        surface_smoothing=2,
        taubin_iter=100,
        dilate=0,
    )
    pv.save_meshio(f"{surface_dir}/ventricles.stl", ventricles)


def process_cerebral_surfaces(
    fs_dir: Path,
    surface_dir: Path,
):
    logger.info(f"Processing surfaces from {fs_dir} -> {surface_dir}")
    Path(surface_dir).mkdir(exist_ok=True)
    seg_mri = sm.load_mri(Path(fs_dir) / "mri/aseg.mgz", dtype=np.int16)
    subcortical_gm = extract_subcortical_gm(seg_mri)
    pv.save_meshio(f"{surface_dir}/subcortical_gm.stl", subcortical_gm)

    for surf in ["lh_pial.stl", "rh_pial.stl", "lh_white.stl", "rh_white.stl"]:
        if not Path(f"{surface_dir}/{surf}").exists():
            try:
                fs_surf_to_stl(
                    f"{fs_dir}/surf", Path(surface_dir), suffix="", verbose=True
                )
            except Exception as e:
                raise RuntimeError(
                    f"Couldn't find {surf} in {surface_dir}. Attempt at creation failed with error:\n{e}"
                )
    preprocess_white_matter_surfaces(fs_dir, surface_dir)
    preprocess_pial_surfaces(surface_dir)


def generate_mesh(
    surface_dir: Path, output: Path, resolution: float, keep_ventricles: bool
):
    logger.info(
        f"Creating cerebral mesh from {surface_dir} -> {output} with resolution {resolution}"
    )
    surface_names = [
        "rh_pial_refined",
        "lh_pial_refined",
        "subcortical_gm",
        "white",
        "ventricles",
    ]
    surface_files = {surf: surface_dir / f"{surf}.stl" for surf in surface_names}
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
    surfaces = [svmtk_surfaces[surf] for surf in surface_names]

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
    pytime.sleep(10)  # Wait another 10 seconds in case of IO delay.
    xdmfdir = output.parent / "mesh_xdmfs"
    xdmfdir.mkdir(exist_ok=True, parents=True)
    mesh2xdmf(str(output.with_suffix(".mesh")), xdmfdir, dim=3)
    _ = xdmf2hdf(xdmfdir, output)


def preprocess_pial_surfaces(
    surface_dir,
    max_edge_length=1.0,
    gapsize=0.2,
):
    logger.info("Preprocessing pial surfaces")
    lh_pial_svm = pyvista2svmtk(
        embed_overlapping(
            pv.read(f"{surface_dir}/lh_pial.stl"),
            pv.read(f"{surface_dir}/lh_white.stl"),
            step_size=2.0,
        )
    )

    lh_pial_svm.separate_close_vertices()
    lh_pial_svm.separate_narrow_gaps(-abs(gapsize))
    lh_pial_svm = taubin_svmtk(lh_pial_svm, 100, 0.1)

    rh_pial_svm = pyvista2svmtk(
        embed_overlapping(
            pv.read(f"{surface_dir}/rh_pial.stl"),
            pv.read(f"{surface_dir}/rh_white.stl"),
            step_size=1.0,
        )
    )
    rh_pial_svm.separate_close_vertices()
    rh_pial_svm.separate_narrow_gaps(-abs(gapsize))
    rh_pial_svm = taubin_svmtk(rh_pial_svm, 100, 0.1)

    svmtk.separate_overlapping_surfaces(lh_pial_svm, rh_pial_svm, edge_movement=-1.0)
    svmtk.separate_close_surfaces(lh_pial_svm, rh_pial_svm, edge_movement=-1.0)

    lh_pial_svm.isotropic_remeshing(max_edge_length, 1, True)
    lh_pial_svm = repair_svmtk(lh_pial_svm)

    rh_pial_svm.isotropic_remeshing(max_edge_length, 1, True)
    rh_pial_svm = repair_svmtk(rh_pial_svm)

    lh_pial_svm.save(f"{surface_dir}/lh_pial_refined.stl")
    rh_pial_svm.save(f"{surface_dir}/rh_pial_refined.stl")


def preprocess_white_matter_surfaces(
    fs_dir,
    surface_dir,
    max_edge_length=0.5,
    remesh_iter=1,
    taubin_iter=100,
    taubin_pass_band=0.1,
):
    logger.info("Preparing white matter surfaces.")
    seg_mri = sm.load_mri(Path(fs_dir) / "mri/aseg.mgz", dtype=np.int16)

    connective = repair_triangulation(
        grow_white_connective_tissue(seg_mri, cc_radius=3, lv_radius=2, smoothing=4)
    )
    svm_surfaces = [
        pyvista2svmtk(repair_triangulation(pv.read(f"{surface_dir}/lh_white.stl"))),
        pyvista2svmtk(connective),
        pyvista2svmtk(repair_triangulation(pv.read(f"{surface_dir}/rh_white.stl"))),
    ]
    white_svm = functools.reduce(surface_union, svm_surfaces)
    white_svm = pyvista2svmtk(repair_triangulation(svmtk2pyvista(white_svm)))
    white_svm.isotropic_remeshing(max_edge_length, remesh_iter, True)
    white = svmtk2pyvista(white_svm)
    white.smooth_taubin(
        taubin_iter, taubin_pass_band, normalize_coordinates=True, inplace=True
    )
    pv.save_meshio(f"{surface_dir}/white.stl", repair_triangulation(white))


def extract_subcortical_gm(
    seg_mri: sm.SimpleMRI, label_smoothing: float = 1.0, binary_smoothing: float = 1.0
):
    logger.info("Extracting subcortical gray matter surfaces.")
    subcort_gm_mask = segmentation_smoothing(
        seg_mri.data, label_smoothing, cutoff_score=1 / 2
    )
    subcortical_gm_mask = np.isin(
        subcort_gm_mask["labels"], default_segmentation_groups()["subcortical-gm"]
    )
    subcortical_gm = binary_image_surface_extraction(
        skimage.morphology.binary_dilation(
            subcortical_gm_mask, skimage.morphology.cube(3)
        ),
        sigma=binary_smoothing,
        cutoff=1 / 2,
    ).transform(seg_mri.affine)
    subcortical_svm = pyvista2svmtk(subcortical_gm)
    subcortical_gm = repair_disconnected_triangulation(svmtk2pyvista(subcortical_svm))
    return subcortical_gm


if __name__ == "__main__":
    meshgen()
