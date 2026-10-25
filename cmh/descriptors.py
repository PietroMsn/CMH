"""Spectral descriptors, port of UtilsFM/waveKernel*.m and UtilsFM/Fmap/hks.m."""
import numpy as np


def _wks_filters(evals, num_times):
    log_e = np.log(np.abs(evals))
    emin, emax = log_e[1], log_e[-1]
    s = 7 * (emax - emin) / num_times
    es = np.linspace(emin + 2 * s, emax - 2 * s, num_times)
    return np.exp(-(log_e[:, None] - es[None, :]) ** 2 / (2 * s ** 2))


def wave_kernel_signature(evecs, evals, A, num_times):
    T = _wks_filters(evals, num_times)
    D = evecs.T @ (A @ evecs ** 2)
    return evecs @ (D @ T)


def wave_kernel_map(evecs, evals, A, num_times, landmarks):
    """WKM centered at each landmark, concatenated along the columns."""
    T = _wks_filters(evals, num_times)
    # A is not used, as in the MATLAB code: the landmark delta is projected
    # with the plain transpose of the basis.
    return np.column_stack([evecs @ (T * evecs[lm][:, None]) for lm in landmarks])


def heat_kernel_signature(evecs, evals, A, ts):
    D = evecs.T @ (A @ evecs ** 2)
    return evecs @ (D @ np.exp(-np.abs(evals)[:, None] * np.asarray(ts)[None, :]))
