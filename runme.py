"""CMH: Coordinates Manifold Harmonics for Functional Remeshing (3DOR 2019).

Python port of runme.m: remeshes the shape in `path_to_remesh` with the
connectivity of the shape in `path_connectivity`. Landmarks are 0-based
(the MATLAB indices minus 1).
"""
import argparse
import time
from pathlib import Path

import igl
import matplotlib.pyplot as plt
import numpy as np

from cmh import (DescriptorOptions, compute_cmh, compute_map_cmh, get_shape,
                 perform_functional_remeshing, scaling_shapes)
from cmh.plotting import plot_cmh, plot_comparison_mh_cmh, plot_results

ROOT = Path(__file__).resolve().parent

path_to_remesh = ROOT / 'shapes' / 'TposeManBad'
path_connectivity = ROOT / 'shapes' / 'TposeWomanFat'

landmarks1 = np.array([1581, 851, 741, 2362, 2782]) - 1
landmarks2 = np.array([5638, 2336, 1066, 8266, 8599]) - 1
landmarks = np.column_stack([landmarks1, landmarks2])

# parameters
cpd_on = False          # refine the map on the hands with Coherent Point Drift
kM1, kM2 = 50, 50       # number of eigenfunctions used (excluding CMH)
desc = DescriptorOptions(wks=True, hks=False, land=True, skip_size=10)
hands_radius = 0.13
a_arap = 0.2
delta_t = 0.03


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--no-plot', action='store_true', help='do not show the plots')
    parser.add_argument('--cpd', action='store_true', help='enable the CPD module')
    parser.add_argument('--out', default=ROOT / 'results' / 'remeshed.obj', type=Path,
                        help='where to save the remeshed shape')
    args = parser.parse_args()
    visualize = not args.no_plot
    times = {}

    M1, L1 = get_shape(path_to_remesh, landmarks1)
    M2, L2 = get_shape(path_connectivity, landmarks2)

    # scaling the two input meshes with the area of the SMPL model
    M1, M2 = scaling_shapes(M1, M2)

    print('Computing CMH for M1 and M2...')
    t = time.perf_counter()
    M1_cmh = compute_cmh(M1, L1.evecs[:, :kM1])
    M2_cmh = compute_cmh(M2, L2.evecs[:, :kM2])
    times['cmh'] = time.perf_counter() - t
    print(f"done in {times['cmh']:f} seconds")

    print('Estimating a functional map in the CMH basis...')
    t = time.perf_counter()
    _, cmh_p2p, _ = compute_map_cmh(M1_cmh, M2_cmh, landmarks, kM1 + 3, kM2 + 3, L1, L2, desc)
    times['map'] = time.perf_counter() - t
    print(f"done in {times['map']:f} seconds")

    if visualize:
        plot_cmh(M1, M2, M1_cmh, M2_cmh, kM1, kM2)
        # the same map with MH only (without the CMH functions)
        _, mh_p2p, _ = compute_map_cmh(M1_cmh, M2_cmh, landmarks, kM1, kM2, L1, L2, desc)
        plot_comparison_mh_cmh(M1, M2, mh_p2p, cmh_p2p)

    _, new2, remesh_times = perform_functional_remeshing(
        M1, M2, landmarks, cmh_p2p, hands_radius, a_arap, delta_t, cpd_on or args.cpd)
    times.update(remesh_times)
    print(f'Total execution time: {sum(times.values()):f} seconds')

    args.out.parent.mkdir(parents=True, exist_ok=True)
    igl.write_triangle_mesh(str(args.out), new2, M2.triv)
    print(f'Remeshed shape saved in {args.out}')

    if visualize:
        plot_results(M1, M2, new2)
        plt.show()


if __name__ == '__main__':
    main()
