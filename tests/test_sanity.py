"""Checks of the Python port that do not need MATLAB."""
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from cmh import compute_cmh, get_shape
from cmh.arap import ArapGradient, fit_rotations
from cmh.fmap import DescriptorOptions, compute_map_cmh

ROOT = Path(__file__).resolve().parent.parent
LANDMARKS = np.column_stack([np.array([1581, 851, 741, 2362, 2782]) - 1,
                             np.array([5638, 2336, 1066, 8266, 8599]) - 1])


@pytest.fixture(scope='module')
def shapes():
    M1, L1 = get_shape(ROOT / 'shapes' / 'TposeManBad', LANDMARKS[:, 0], k=60)
    M2, L2 = get_shape(ROOT / 'shapes' / 'TposeWomanFat', LANDMARKS[:, 1], k=60)
    return M1, L1, M2, L2


def test_lb_basis_is_orthonormal(shapes):
    _, L1, _, _ = shapes
    gram = L1.evecs.T @ (L1.A @ L1.evecs)
    assert np.allclose(gram, np.eye(gram.shape[0]), atol=1e-8)
    assert np.all(np.diff(L1.evals) >= 0) and abs(L1.evals[0]) < 1e-6


def test_cmh_is_orthonormal_and_spans_coordinates(shapes):
    M1, L1, _, _ = shapes
    cmh = compute_cmh(M1, L1.evecs[:, :50])
    gram = cmh.T @ (M1.Ae @ cmh)
    assert np.allclose(gram, np.eye(53), atol=1e-8)
    # the coordinates are exactly reconstructed by the CMH basis
    rec = cmh @ (cmh.T @ (M1.Ae @ M1.vert))
    assert np.allclose(rec, M1.vert, atol=1e-8)


def test_arap_gradient_vanishes_at_rotated_rest_pose(shapes):
    M1, _, _, _ = shapes
    grad = ArapGradient(M1.vert, M1.triv)
    Q = Rotation.from_euler('xyz', [0.3, -0.7, 1.1]).as_matrix()
    G = grad(M1.vert @ Q.T + [1.0, 2.0, 3.0])
    assert np.abs(G).max() < 1e-8 * np.abs(grad.L @ M1.vert).max() + 1e-10


def test_arap_gradient_matches_finite_differences(shapes):
    M1, _, _, _ = shapes
    grad = ArapGradient(M1.vert, M1.triv)
    V = M1.vert
    rng = np.random.default_rng(0)
    U = V + 0.01 * rng.standard_normal(V.shape)

    def energy(U):
        # E as in arap_gradient.m, with the best-fit rotations of U
        R = fit_rotations(np.stack([kt @ U for kt in grad.KT], axis=1))
        dV = sum(grad.K[c] @ R[:, :, c] for c in range(3))
        return np.trace(-U.T @ (0.5 * grad.L @ U) - U.T @ dV - V.T @ (0.5 * grad.L @ V))

    D = rng.standard_normal(V.shape)
    eps = 1e-6
    fd = (energy(U + eps * D) - energy(U - eps * D)) / (2 * eps)
    assert fd == pytest.approx(np.sum(grad(U) * D), rel=1e-4)


def test_cmh_map_is_accurate_on_landmarks(shapes):
    M1, L1, M2, L2 = shapes
    cmh1 = compute_cmh(M1, L1.evecs[:, :50])
    cmh2 = compute_cmh(M2, L2.evecs[:, :50])
    _, p2p, _ = compute_map_cmh(cmh1, cmh2, LANDMARKS, 53, 53, L1, L2, DescriptorOptions())
    assert p2p.shape == (len(M2.vert),)
    # the landmarks of M2 land close to the ones of M1
    err = np.linalg.norm(M1.vert[p2p[LANDMARKS[:, 1]]] - M1.vert[LANDMARKS[:, 0]], axis=1)
    diameter = np.ptp(M1.vert, axis=0).max()
    assert np.all(err < 0.05 * diameter)
