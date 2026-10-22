import numpy as np
import pytest

from dti.reslice_dti import eigenvector_rotation


def rotation(axis: int, degrees: float) -> np.ndarray:
    a = np.radians(degrees)
    i, j = [k for k in range(3) if k != axis]
    R = np.eye(3)
    R[i, i], R[i, j], R[j, i], R[j, j] = np.cos(a), -np.sin(a), np.sin(a), np.cos(a)
    return R


def affine(directions: np.ndarray, zooms=(2.5, 2.0, 3.0)) -> np.ndarray:
    A = np.eye(4)
    A[:3, :3] = directions @ np.diag(zooms)
    A[:3, 3] = [100.0, -80.0, -40.0]
    return A


def fsl_vector(world_direction: np.ndarray, A: np.ndarray) -> np.ndarray:
    """Unit vector as dtifit writes it: voxel-index differences scaled by the voxel size,
    with the first axis reversed for a positive-determinant affine (FSL's radiological convention).
    """
    p = np.array([10.0, 20.0, 30.0])
    Ainv = np.linalg.inv(A)
    ijk = [Ainv[:3, :3] @ q + Ainv[:3, 3] for q in (p, p + world_direction)]
    v = (ijk[1] - ijk[0]) * np.linalg.norm(A[:3, :3], axis=0)
    if np.linalg.det(A[:3, :3]) > 0:
        v[0] = -v[0]
    return v / np.linalg.norm(v)


OBLIQUE = rotation(2, 9.0) @ rotation(0, 8.0)
LAS = affine(OBLIQUE @ np.diag([-1.0, 1.0, 1.0]))  # like the Gonzo DWI
RAS = affine(OBLIQUE)
TRANSFORM = np.eye(4)  # fixed world -> moving world, as written by greedy
TRANSFORM[:3, :3] = rotation(0, 5.0) @ rotation(1, -3.0)
TRANSFORM[:3, 3] = [2.0, -0.3, -2.9]
DIRECTIONS = [
    np.array(d, float) / np.linalg.norm(d)
    for d in ([1, 1, 0], [1, 0, -1], [0.2, -1, 0.5])
]


@pytest.mark.parametrize("A", [LAS, RAS], ids=["LAS", "RAS"])
@pytest.mark.parametrize("u", DIRECTIONS)
def test_eigenvectors_to_moving_world(A, u):
    assert np.allclose(eigenvector_rotation(A) @ fsl_vector(u, A), u)


@pytest.mark.parametrize("A", [LAS, RAS], ids=["LAS", "RAS"])
@pytest.mark.parametrize("u_fixed", DIRECTIONS)
def test_eigenvectors_to_fixed_world(A, u_fixed):
    u_moving = TRANSFORM[:3, :3] @ u_fixed
    assert np.allclose(
        eigenvector_rotation(A, TRANSFORM) @ fsl_vector(u_moving, A), u_fixed
    )


def test_rotation_is_orthogonal():
    R = eigenvector_rotation(LAS, TRANSFORM)
    assert np.allclose(R @ R.T, np.eye(3))
