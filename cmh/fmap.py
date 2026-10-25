"""Functional map estimation, port of utils_CMH/compute_map_CMH.m and
UtilsFM/icp_refine.m.

minConf_PQN is replaced by L-BFGS on the free entries of C: the MATLAB
projection only pins the first column, and the energy is a convex quadratic,
so both converge to the same minimizer.
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

from .descriptors import heat_kernel_signature, wave_kernel_map, wave_kernel_signature


@dataclass
class DescriptorOptions:
    wks: bool = True
    hks: bool = False
    land: bool = True
    skip_size: int = 10


def _descriptors(L, landmarks, desc):
    fct = []
    if desc.wks:
        fct.append(wave_kernel_signature(L.evecs, L.evals, L.A, 200))
    if desc.land:
        fct.append(wave_kernel_map(L.evecs, L.evals, L.A, 200, landmarks))
    if desc.hks:
        ts = np.exp(np.linspace(np.log(0.03), np.log(0.5), 100))
        fct.append(np.log(heat_kernel_signature(L.evecs, L.evals, L.A, ts)))
    return np.column_stack(fct)[:, ::desc.skip_size]


def _multiplication_operators(phi, a, fct):
    """Phi' A diag(f) Phi for every descriptor f, stacked as (nf, k, k)."""
    phi_a = phi * a[:, None]
    return np.stack([phi_a.T @ (phi * f[:, None]) for f in fct.T])


def icp_refine(phi1, phi2, C, n_iter):
    k2, k1 = C.shape
    for _ in range(n_iter):
        matches = cKDTree(phi1 @ C.T).query(phi2)[1]
        W = np.linalg.lstsq(phi2, phi1[matches], rcond=None)[0]
        U, _, Vt = np.linalg.svd(W)
        C = U @ np.eye(k2, k1) @ Vt
    return C, matches


def compute_map_cmh(basis1, basis2, landmarks, k1, k2, L1, L2, desc):
    """Functional map from shape 1 to shape 2 and the induced p2p map.

    `landmarks` is (n_landmarks, 2) with 0-based indices. Returns the (k2, k1)
    map after ICP, the vertex map p2p (for every vertex of shape 2, the index
    of its match on shape 1) and the map before ICP.
    """
    phi1, phi2 = basis1[:, :k1], basis2[:, :k2]
    a1, a2 = L1.A.diagonal(), L2.A.diagonal()

    fct_src = _descriptors(L1, landmarks[:, 0], desc)
    fct_tar = _descriptors(L2, landmarks[:, 1], desc)

    # both sides are normalized with the norms on the source
    no = np.sqrt(np.einsum('ij,i,ij->j', fct_src, a1, fct_src))
    fct_src, fct_tar = fct_src / no, fct_tar / no

    op_src = _multiplication_operators(phi1, a1, fct_src)
    op_tar = _multiplication_operators(phi2, a2, fct_tar)
    f_src = phi1.T @ (a1[:, None] * fct_src)
    f_tar = phi2.T @ (a2[:, None] * fct_tar)

    dlb = (L1.evals[:k1, None] - L2.evals[None, :k2]) ** 2
    dlb = dlb / np.linalg.norm(dlb) ** 2
    # MATLAB multiplies vec(C) by vec(Dlb) element-wise even though C is
    # k2 x k1 and Dlb is k1 x k2: keep the same pairing of the entries.
    dlb = dlb.ravel(order='F').reshape((k2, k1), order='F')

    const = np.sign(phi1[0, 0] * phi2[0, 0]) * np.sqrt(a2.sum() / a1.sum())

    a, b, c = 1e-1, 1.0, 1e-3   # descriptors, commutativity with descriptors, with LB

    def unpack(x):
        C = np.zeros((k2, k1))
        C[0, 0] = const
        C[:, 1:] = x.reshape((k2, k1 - 1), order='F')
        return C

    def energy(x):
        C = unpack(x)
        r_desc = C @ f_src - f_tar
        r_comm = op_tar @ C - C @ op_src
        e = (a * np.sum(r_desc ** 2) + b * np.sum(r_comm ** 2) + c * np.sum(C ** 2 * dlb)) / 2
        g = (a * r_desc @ f_src.T
             + b * np.sum(op_tar.transpose(0, 2, 1) @ r_comm - r_comm @ op_src.transpose(0, 2, 1), axis=0)
             + c * C * dlb)
        return e, g[:, 1:].ravel(order='F')

    res = minimize(energy, np.zeros(k2 * (k1 - 1)), jac=True, method='L-BFGS-B',
                   options={'maxiter': 10000, 'maxcor': 50, 'ftol': 1e-15, 'gtol': 1e-10})
    C0 = unpack(res.x)

    C, _ = icp_refine(phi1, phi2, C0, 5)
    p2p = cKDTree(phi1 @ C.T).query(phi2)[1]
    return C, p2p, C0
