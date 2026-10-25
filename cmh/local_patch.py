"""Refinement of the map on local patches with Coherent Point Drift, port of
utils/LocalPat.m and utils/LocalPatch/LocalPatchs/*.

CPD comes from pycpd instead of the CPD2 mex files, and geodesic distances
from libigl's exact geodesics instead of fastmarchmex.
"""
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from itertools import repeat

import igl
import numpy as np
from pycpd import DeformableRegistration, RigidRegistration
from scipy.spatial import cKDTree


@dataclass
class Patch:
    vert: np.ndarray
    triv: np.ndarray
    fullshape_idx: np.ndarray   # index of each patch vertex in the full shape


def extract_patch(vert, triv, vert_mask):
    """Patch made of the triangles touching at least one vertex in the mask."""
    tris = triv[vert_mask[triv].any(axis=1)]
    keep = np.unique(tris)
    old_to_new = np.full(len(vert), -1)
    old_to_new[keep] = np.arange(len(keep))
    return Patch(vert=vert[keep], triv=old_to_new[tris], fullshape_idx=keep)


def geodesic_distances(vert, triv, source):
    return igl.exact_geodesic(vert, triv, np.array([source]), np.array([], dtype=np.int64),
                              np.arange(len(vert)), np.array([], dtype=np.int64))


def geodesic_matrix(vert, triv, samples):
    D = np.stack([geodesic_distances(vert, triv, s)[samples] for s in samples], axis=1)
    return (D + D.T) / 2


def fps_euclidean(vert, n, seed):
    samples = [seed]
    d = np.linalg.norm(vert - vert[seed], axis=1)
    for _ in range(n - 1):
        samples.append(int(np.argmax(d)))
        d = np.minimum(d, np.linalg.norm(vert - vert[samples[-1]], axis=1))
    return np.array(samples)


def _normalize(x):
    mean = x.mean(axis=0)
    scale = np.sqrt(np.sum((x - mean) ** 2) / len(x))
    return (x - mean) / scale, mean, scale


def align_icp_cpd(X, Y, lam=10):
    """Rigid then non-rigid CPD of Y onto X, as align_icp_cpd_hand.m.

    Both point sets are normalized to zero mean and unit variance, as
    cpd_register does, and the result is brought back to the frame of X.
    """
    Xn, x_mean, x_scale = _normalize(X)
    Yn, _, _ = _normalize(Y)
    Yn, _ = RigidRegistration(X=Xn, Y=Yn, w=0.6, max_iterations=200, tolerance=1e-8).register()
    Yn, _, _ = _normalize(Yn * x_scale + x_mean)
    Yn, _ = DeformableRegistration(X=Xn, Y=Yn, w=0.6, alpha=lam, beta=4,
                                   max_iterations=200, tolerance=1e-8).register()
    return Yn * x_scale + x_mean


def local_patch_refine(M1, M2, landmarks, p2p, dist_thresh, n_angles=16):
    """Re-match the vertices of M2 around landmarks[1] to the patch of M1
    around landmarks[0], trying n_angles initial rotations around z."""
    P1 = extract_patch(M1.vert, M1.triv, geodesic_distances(M1.vert, M1.triv, landmarks[0]) < dist_thresh)
    P2 = extract_patch(M2.vert, M2.triv, geodesic_distances(M2.vert, M2.triv, landmarks[1]) < dist_thresh)
    tree = cKDTree(P1.vert)

    n_pts = fps_euclidean(P2.vert, 100, 0)
    DN = geodesic_matrix(P2.vert, P2.triv, n_pts)

    rotated = []
    for angle in np.linspace(0, 2 * np.pi * (1 - 1 / n_angles), n_angles):
        R = np.array([[np.cos(angle), -np.sin(angle), 0],
                      [np.sin(angle), np.cos(angle), 0],
                      [0, 0, 1]])
        rotated.append(P2.vert @ R)
    # the registrations are independent, as in the parfor of LocalPat.m;
    # one BLAS thread per worker, read by the workers when they import numpy
    blas_vars = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')
    saved = {k: os.environ.get(k) for k in blas_vars}
    os.environ.update({k: '1' for k in blas_vars})
    try:
        with ProcessPoolExecutor() as pool:
            solutions = list(pool.map(align_icp_cpd, repeat(P1.vert), rotated, repeat(3)))
    finally:
        for k, v in saved.items():
            if v is None:
                os.environ.pop(k)
            else:
                os.environ[k] = v

    best_matches, best_distortion = None, np.inf
    for deformed in solutions:
        matches = tree.query(deformed)[1]
        DM = geodesic_matrix(P1.vert, P1.triv, matches[n_pts])
        distortion = np.abs(DM - DN).sum()
        if distortion < best_distortion:
            best_matches, best_distortion = matches, distortion

    p2p = p2p.copy()
    p2p[P2.fullshape_idx] = P1.fullshape_idx[best_matches]
    return p2p
