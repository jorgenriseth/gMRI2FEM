import subprocess
import tempfile
from pathlib import Path
from typing import Literal, Optional
from scipy.spatial import KDTree

import pyvista as pv
import SVMTK as svmtk
import dolfin as df
import simple_mri as sm
import skimage
import numpy as np
import pymeshfix
from loguru import logger
from gmri2fem.utils import grow_restricted, largest_island

WHITE = [2, 41]
CEREBRAL_CORTEX = [3, 42]
LATERAL_VENTRICLES = [4, 43]
GENERIC_CSF = [24]
CC = [251, 252, 253, 254, 255]
DC = [28, 60]


def fs_surf_to_stl(
    fs_surface_dir: Path | str,
    output_dir: Path,
    suffix: Optional[str] = None,
    verbose: bool = False,
):
    suffix = "" if suffix is None else suffix

    fs_surfaces = ["rh.pial", "lh.pial", "rh.white", "lh.white"]
    logger.info("Converting FS-surfaces to .stl")
    for surface in fs_surfaces:
        input = f"{fs_surface_dir}/{surface}"
        output = (output_dir / (surface.replace(".", "_") + suffix)).with_suffix(".stl")
        redirect = ">> /dev/null" if not verbose else ""
        subprocess.run(
            f"mris_convert --to-scanner {input} {output} {redirect}",
            shell=True,
        ).check_returncode()


def print_edge_info(mesh: pv.UnstructuredGrid):
    edges = mesh.extract_all_edges()
    edge_lengths = edges.compute_cell_sizes(length=True)["Length"]  # type: ignore
    print(f"Min edge length: {edge_lengths.min()}")
    print(f"Max edge length: {edge_lengths.max()}")
    print(f"Mean edge length: {edge_lengths.mean()}")


def replace_at_index(s: str, idx: int, value: Literal["0", "1"]) -> str:
    return s[:idx] + value + s[idx + 1 :]


def expand_subdomain_string(smap_string: str) -> list[str]:
    idx = smap_string.find(".")
    if idx == -1:
        return [smap_string]
    else:
        return expand_subdomain_string(
            replace_at_index(smap_string, idx, "0")
        ) + expand_subdomain_string(replace_at_index(smap_string, idx, "1"))


def subdomain_mapper(smap: svmtk.SubdomainMap, smap_string: str, tag: int):
    for s in expand_subdomain_string(smap_string):
        smap.add(s, tag)
    return smap


def map_subdomains_to_boundaries(mesh, subdomains):
    mesh.init()
    d = mesh.topology().dim()
    boundaries = df.MeshFunction("size_t", mesh, d - 1, 2)
    facet_labels = np.zeros_like(boundaries.array())
    for facet in df.facets(mesh):
        ent = facet.entities(d)
        if len(ent) == 1:
            facet_labels[facet.index()] = subdomains.array()[ent[0]]
    boundaries.array()[:] = facet_labels
    return boundaries


def map_boundaries_to_subdomains(mesh, boundaries):
    mesh.init()
    d = mesh.topology().dim()
    subdomains = df.MeshFunction("size_t", mesh, d, 0)
    cell_labels = np.zeros_like(subdomains.array())
    boundary_facets = [
        facet.index() for facet in df.facets(mesh) if len(facet.entities(d)) == 1
    ]
    for cell in df.cells(mesh):
        cell_facets = cell.entities(d - 1)
        cell_boundary_facets = cell_facets[np.isin(cell_facets, boundary_facets)]
        if len(cell_boundary_facets) > 0:
            cell_boundary_tags = boundaries.array()[:][cell_boundary_facets]
            values, counts = np.unique(cell_boundary_tags, return_counts=True)
            cell_labels[cell.index()] = values[counts.argmax()]
    subdomains.array()[:] = cell_labels
    return subdomains


def tagged_facet_extraction(input: Path, output):
    # Load the tetrahedral mesh
    mesh = pv.read(input)  # Or create mesh here
    if mesh is None:
        raise ValueError(f"Can't load mesh from {input}")

    # Assuming cell tags are stored in 'cell_tags'
    cell_tags = mesh.cell_data["subdomains"]

    # Check if all cells are tetrahedra
    from vtk import VTK_TETRA

    if not np.all(mesh.celltypes == VTK_TETRA):
        raise ValueError("The mesh must contain only tetrahedral cells.")

    # Extract cell connectivity
    cells = mesh.cells.reshape(-1, 5)
    assert np.all(cells[:, 0] == 4)  # Ensure all cells have 4 points
    cell_point_ids = cells[:, 1:]

    # Build face-cell mapping
    face_dict = {}
    for cell_id, cell_points in enumerate(cell_point_ids):
        faces = [
            tuple(sorted([cell_points[0], cell_points[1], cell_points[2]])),
            tuple(sorted([cell_points[0], cell_points[1], cell_points[3]])),
            tuple(sorted([cell_points[0], cell_points[2], cell_points[3]])),
            tuple(sorted([cell_points[1], cell_points[2], cell_points[3]])),
        ]
        for face in faces:
            if face in face_dict:
                face_dict[face]["cells"].append(cell_id)
            else:
                face_dict[face] = {"cells": [cell_id]}

    # Assign tags to faces
    subdomain_pair_to_id = {}
    current_id = cell_tags.max() + 1  # Start from 1 since 0 is reserved
    face_data = []

    for face, data in face_dict.items():
        cells_adjacent = data["cells"]
        if len(cells_adjacent) == 1:
            # Boundary face
            cell_id = cells_adjacent[0]
            face_tag = cell_tags[cell_id]
        elif len(cells_adjacent) == 2:
            # Internal face
            cell_id1, cell_id2 = cells_adjacent
            tag1 = cell_tags[cell_id1]
            tag2 = cell_tags[cell_id2]
            if tag1 == tag2:
                face_tag = 0
            else:
                subdomain_pair = tuple(sorted((tag1, tag2)))
                if subdomain_pair not in subdomain_pair_to_id:
                    subdomain_pair_to_id[subdomain_pair] = current_id
                    current_id += 1
                face_tag = subdomain_pair_to_id[subdomain_pair]
        else:
            raise ValueError(
                "A face is adjacent to more than two cells, which should not happen."
            )
        face_data.append({"face": face, "tag": face_tag})

    # Prepare faces and tags for PolyData
    faces_list = []
    face_tags = []
    for face_info in face_data:
        face = face_info["face"]
        face_tag = face_info["tag"]
        faces_list.extend([3, face[0], face[1], face[2]])
        face_tags.append(face_tag)

    faces_array = np.array(faces_list)
    face_tags = np.array(face_tags)

    # Create the face mesh
    face_mesh = pv.PolyData(mesh.points, faces_array)
    face_mesh.cell_data["face_tags"] = face_tags

    # Save the result
    face_mesh.save(output)


def repair_triangulation(mesh):
    if isinstance(mesh, pv.PolyData):
        v, f = mesh.points, mesh.faces.reshape((-1, 4))[:, 1:]
    else:
        v, f = mesh.points, mesh.cells.reshape((-1, 4))[:, 1:]
    meshfix = pymeshfix.MeshFix(v, f)
    meshfix.repair()

    num_points_stack = 3 * np.ones(
        (meshfix.faces.shape[0], 1), dtype=meshfix.faces.dtype
    )
    pv_faces = np.concatenate((num_points_stack, meshfix.faces), axis=1).ravel()
    return pv.PolyData(meshfix.v, pv_faces)


def repair_disconnected_triangulation(surfaces):
    surface_components = surfaces.connectivity()
    region_ids = surface_components.cell_data["RegionId"]
    unique_regions = np.unique(region_ids)
    region_meshes = []
    for rid in unique_regions:
        region_mask = region_ids == rid
        cell_indices = np.where(region_mask)[0]
        region_mesh = surface_components.extract_cells(cell_indices)
        region_meshes.append(repair_triangulation(region_mesh.triangulate()))
    return region_meshes[0].merge(region_meshes[1:])


def surface_union(x, y):
    x.union(y)
    return x


def taubin_smooth(svm_surface, taubin_iter, taubin_pass_band):
    pv_surface = svmtk2pyvista(svm_surface)
    pv_surface.smooth_taubin(
        taubin_iter, taubin_pass_band, normalize_coordinates=True, inplace=True
    )
    return pyvista2svmtk(pv_surface)


def pyvista2svmtk(
    pv_grid: pv.DataObject, suffix: Optional[str] = None
) -> svmtk.Surface:
    ft = ".stl" if suffix is None else suffix
    with tempfile.TemporaryDirectory() as tmp_path:
        tmpfile = Path(tmp_path) / f"tmpsurf{ft}"
        pv.save_meshio(tmpfile, pv_grid)
        svmtk_grid = svmtk.Surface(str(tmpfile))
    return svmtk_grid


def svmtk2pyvista(
    svmtk_surface: svmtk.Surface, suffix: Optional[str] = None
) -> pv.DataObject:
    ft = ".stl" if suffix is None else suffix
    with tempfile.TemporaryDirectory() as tmp_path:
        tmpfile = Path(tmp_path) / f"tmpsurf{ft}"
        svmtk_surface.save(str(tmpfile))
        pv_grid = pv.read(tmpfile)
    return pv_grid


def repair_svmtk(surface):
    return pyvista2svmtk(repair_triangulation(svmtk2pyvista(surface)))


def taubin_svmtk(surface, n_iter, pass_band):
    pv_surface = svmtk2pyvista(surface)
    pv_surface.smooth_taubin(
        n_iter, pass_band, normalize_coordinates=True, inplace=True
    )
    return pyvista2svmtk(pv_surface)


def embed_overlapping(pial, white, step_size=1.0):
    tree = KDTree(white.points)
    d_kdtree, _ = tree.query(pial.points)
    nonzero_distance_verts = np.argwhere(d_kdtree != 0.0)
    warped = pial.compute_normals(point_normals=True).warp_by_vector(
        "Normals", factor=-step_size
    )
    warped.points[nonzero_distance_verts] = pial.points[nonzero_distance_verts]
    return repair_triangulation(warped)


def remove_ventricles(
    mesh_path: Path, output_path: Path, ventricles_label: Optional[int] = None
):
    mesh = pv.read_meshio(mesh_path)
    all_labels = np.unique(mesh.cell_data["label"])
    if ventricles_label is None:
        ventricles_label = all_labels.max()
    brain_labels = all_labels[all_labels != ventricles_label]
    brain = mesh.extract_values(brain_labels)
    pv.save_meshio(output_path, brain)


def grow_white_connective_tissue(
    seg_mri: sm.SimpleMRI, cc_radius: int = 2, lv_radius: int = 3, smoothing: float = 2
):
    CC_mask = grow_restricted(
        np.isin(seg_mri.data, CC), np.isin(seg_mri.data, WHITE + [24]), cc_radius
    )
    LV_mask = grow_restricted(
        np.isin(seg_mri.data, LATERAL_VENTRICLES),
        ~np.isin(seg_mri.data, CEREBRAL_CORTEX),
        lv_radius,
    )
    surf = binary_image_surface_extraction(
        largest_island(LV_mask) + CC_mask, sigma=smoothing
    )
    return surf.transform(seg_mri.affine)


def binary_image_surface_extraction(
    vol: np.ndarray,
    sigma: float = 0,
    cutoff=0.5,
    to_cell: bool = False,
    connected: bool = False,
) -> pv.PolyData:
    grid = pv.ImageData(dimensions=vol.shape, spacing=[1] * 3, origin=[0] * 3)
    if sigma > 0:
        grid["labels"] = skimage.filters.gaussian(vol, sigma=sigma).ravel(order="F")
    else:
        grid["labels"] = vol.ravel(order="F")
    if to_cell:
        thresh = grid.points_to_cells().threshold(cutoff)
        surf = thresh.extract_surface()  # type: ignore
    else:
        surf = grid.contour([cutoff])  # type: ignore
    surf.clear_data()  # type: ignore
    return surf
