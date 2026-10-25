"""Functional remeshing, port of perform_functional_remeshing.m and
utils/arap_dato.m."""
import time

import igl
import numpy as np

from .arap import ArapGradient
from .local_patch import local_patch_refine


def _closest_points(vert, triv, query):
    return igl.point_mesh_squared_distance(query, vert, triv)[2]


def arap_dato(new, M1, M2, delta_t, a_arap, n_iter=399):
    """ARAP flow towards M2 plus attraction to the surface of M1."""
    arap_grad = ArapGradient(M2.vert, M2.triv)
    surface_points = _closest_points(M1.vert, M1.triv, new)
    for i in range(1, n_iter + 1):
        if i % 50 == 0:
            print(f'iter:{i}')
            surface_points = _closest_points(M1.vert, M1.triv, new)
        new = new - delta_t * (a_arap * arap_grad(new) + (new - surface_points))
    return new


def perform_functional_remeshing(M1, M2, landmarks, p2p, hands_radius, a_arap,
                                 delta_t, cpd_on=False):
    """Geometry of M1 with the connectivity of M2.

    Returns the vertices after the ARAP flow, the final vertices (both to use
    with M2.triv) and the time spent in each step.
    """
    times = {}

    t = time.perf_counter()
    if cpd_on:
        print('Applying Coherent Point Drift on local patches...')
        p2p = local_patch_refine(M1, M2, landmarks[1], p2p, hands_radius)
        p2p = local_patch_refine(M1, M2, landmarks[4], p2p, hands_radius)
    new = M1.vert[p2p]
    times['cpd'] = time.perf_counter() - t
    print(f"done in {times['cpd']:f} seconds")

    print('processing arap for the mesh consistency...')
    t = time.perf_counter()
    arap_grad = ArapGradient(M2.vert, M2.triv)
    for _ in range(200):
        new = new - delta_t * a_arap * arap_grad(new)
    times['arap'] = time.perf_counter() - t
    print(f"done in {times['arap']:f} seconds")

    print('Minimizing the ARAP energy to fit the original data...')
    t = time.perf_counter()
    new2 = arap_dato(new, M1, M2, delta_t, a_arap)
    times['arap_dato'] = time.perf_counter() - t
    print(f"done in {times['arap_dato']:f} seconds")

    return new, new2, times
