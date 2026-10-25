"""Matplotlib versions of the plots in utils_CMH/."""
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def _draw(ax, vert, triv, f=None, edges=False, title=None):
    polys = vert[triv]
    if f is None:
        colors = np.array([0.8, 0.8, 0.8, 1.0])
    else:
        fv = f[triv].mean(axis=1)
        lim = np.abs(fv).max()
        colors = plt.get_cmap('bwr')((fv / lim + 1) / 2)
    ax.add_collection3d(Poly3DCollection(polys, facecolors=colors,
                                         edgecolors='k' if edges else 'none',
                                         linewidths=0.1, shade=f is None))
    lo, hi = vert.min(axis=0), vert.max(axis=0)
    ax.set(xlim=(lo[0], hi[0]), ylim=(lo[1], hi[1]), zlim=(lo[2], hi[2]))
    ax.set_box_aspect(hi - lo, zoom=1.4)
    ax.view_init(elev=90, azim=-90)   # MATLAB view([0 90])
    ax.set_axis_off()
    if title:
        ax.set_title(title)


def plot_cmh(M1, M2, M1_cmh, M2_cmh, k1, k2):
    fig = plt.figure(figsize=(12, 8))
    for row, (M, cmh, k) in enumerate([(M1, M1_cmh, k1), (M2, M2_cmh, k2)]):
        for i in range(3):
            ax = fig.add_subplot(2, 3, 3 * row + i + 1, projection='3d')
            _draw(ax, M.vert, M.triv, cmh[:, k + i])
    return fig


def plot_comparison_mh_cmh(M1, M2, mh_p2p, cmh_p2p):
    fig = plt.figure(figsize=(15, 6))
    meshes = [(M1.vert, M1.triv, 'Geometry'),
              (M1.vert[mh_p2p], M2.triv, 'MH reconstruction k=50'),
              (M1.vert[cmh_p2p], M2.triv, 'CMH reconstruction k=53')]
    for i, (vert, triv, title) in enumerate(meshes):
        _draw(fig.add_subplot(1, 3, i + 1, projection='3d'), vert, triv, edges=True, title=title)
    return fig


def plot_results(M1, M2, new2):
    fig = plt.figure(figsize=(15, 6))
    meshes = [(M1.vert, M1.triv, 'Geometry'),
              (M2.vert, M2.triv, 'Connectivity'),
              (new2, M2.triv, 'Geometry + connectivity')]
    for i, (vert, triv, title) in enumerate(meshes):
        _draw(fig.add_subplot(1, 3, i + 1, projection='3d'), vert, triv, edges=True, title=title)
    return fig
