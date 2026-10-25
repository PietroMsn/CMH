"""Coordinates Manifold Harmonics, port of compute_CMH.m."""
import numpy as np


def compute_cmh(shape, basis):
    """Extend `basis` with the 3 CMH functions of `shape`.

    Returns a (n, k + 3) matrix: the k input functions followed by the
    Ae-orthonormalized residuals of the x, y, z coordinates.
    """
    Ae = shape.Ae
    error = shape.vert - basis @ (basis.T @ (Ae @ shape.vert))

    cmh = basis
    for i in range(3):
        f = error[:, i] - cmh @ (cmh.T @ (Ae @ error[:, i]))
        f = f / np.sqrt(f @ (Ae @ f))
        cmh = np.column_stack([cmh, f])
    return cmh
