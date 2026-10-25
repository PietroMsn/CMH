"""Step by step comparison with the MATLAB pipeline.

Needs tests/reference/reference_cmh.mat, written by matlab/export_reference.m.
Every step starts from the MATLAB output of the previous one, so that the
sign ambiguity of the eigenvectors does not propagate.
"""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import scipy.io as sio
import scipy.sparse as sp

from cmh.arap import ArapGradient
from cmh.basis import compute_cmh
from cmh.fmap import DescriptorOptions, compute_map_cmh
from cmh.geometry import LBBasis, get_shape, scaling_shapes
from cmh.remeshing import perform_functional_remeshing

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / 'tests' / 'reference' / 'reference_cmh.mat'
LANDMARKS = np.column_stack([np.array([1581, 851, 741, 2362, 2782]) - 1,
                             np.array([5638, 2336, 1066, 8266, 8599]) - 1])

pytestmark = pytest.mark.skipif(not REFERENCE.exists(),
                                reason='run matlab/export_reference.m to create the reference')


@pytest.fixture(scope='module')
def ref():
    r = sio.loadmat(str(REFERENCE), squeeze_me=True)
    for key in ('M1_TRIV', 'M2_TRIV', 'CMH_p2p', 'MH_p2p'):
        r[key] = r[key].astype(np.int64) - 1
    return r


def _shape(ref, i):
    A = sp.diags(ref[f'L{i}_A']).tocsr()
    M = SimpleNamespace(vert=ref[f'M{i}_VERT'], triv=ref[f'M{i}_TRIV'],
                        A=sp.diags(ref[f'M{i}_A']).tocsr(), Ae=A)
    L = LBBasis(evecs=ref[f'L{i}_evecs'], evals=ref[f'L{i}_evals'], A=A, W=None)
    return M, L


def _rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def test_shapes_and_lb_basis(ref):
    M1, L1 = get_shape(ROOT / 'shapes' / 'TposeManBad', LANDMARKS[:, 0])
    M2, L2 = get_shape(ROOT / 'shapes' / 'TposeWomanFat', LANDMARKS[:, 1])
    M1, M2 = scaling_shapes(M1, M2)
    for i, (M, L) in enumerate([(M1, L1), (M2, L2)], start=1):
        assert np.array_equal(M.triv, ref[f'M{i}_TRIV'])
        assert _rel(M.vert, ref[f'M{i}_VERT']) < 1e-10
        assert _rel(L.A.diagonal(), ref[f'L{i}_A']) < 1e-10
        assert _rel(L.evals, ref[f'L{i}_evals']) < 1e-6
        # eigenvectors match up to sign: compare the first 50, used by CMH
        overlap = np.abs(np.einsum('ij,i,ij->j', L.evecs[:, :50], ref[f'L{i}_A'],
                                   ref[f'L{i}_evecs'][:, :50]))
        assert np.all(overlap > 0.99)


def test_cmh(ref):
    for i in (1, 2):
        M, L = _shape(ref, i)
        assert _rel(compute_cmh(M, L.evecs[:, :50]), ref[f'M{i}_CMH']) < 1e-8


@pytest.mark.parametrize('name,k', [('CMH', 53), ('MH', 50)])
def test_functional_map(ref, name, k):
    _, L1 = _shape(ref, 1)
    _, L2 = _shape(ref, 2)
    C, p2p, C0 = compute_map_cmh(ref['M1_CMH'], ref['M2_CMH'], LANDMARKS, k, k, L1, L2,
                                 DescriptorOptions())
    assert _rel(C0, ref[f'{name}_C0']) < 1e-3
    assert np.mean(p2p == ref[f'{name}_p2p']) > 0.95


def test_arap_gradient(ref):
    G = ArapGradient(ref['M2_VERT'], ref['M2_TRIV'])(ref['arap_U'])
    assert _rel(G, ref['arap_G']) < 1e-6


def test_remeshing(ref):
    M1, _ = _shape(ref, 1)
    M2, _ = _shape(ref, 2)
    new, new2, _ = perform_functional_remeshing(M1, M2, LANDMARKS, ref['CMH_p2p'],
                                                0.13, 0.2, 0.03, cpd_on=False)
    assert _rel(new, ref['new_VERT']) < 1e-6
    assert _rel(new2, ref['new2_VERT']) < 1e-4
