from __future__ import with_statement

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os

if __name__ == "__main__":

    fig_dir = "plots"
    data_dir = "z_dist"

    d_z = None
    z_min = None
    z_max = None

    list_of_files = os.listdir(data_dir)
    dtype = np.dtype([('z', float), ('ct', float)])
    for file_name in list_of_files:
        if 'stars' in file_name:
            data = np.genfromtxt(os.path.join(data_dir, file_name), dtype=dtype)
            dd = np.diff(data['z']).min()
            local_min = data['z'].min()
            local_max = data['z'].max()
            if d_z is None or dd<d_z:
                d_z=dd
            if z_min is None or local_min<z_min:
                z_min = local_min
            if z_max is None or local_max>z_max:
                z_max = local_max


    z_grid = np.arange(z_min, z_max+0.5*d_z, d_z)
    ct_grid = np.zeros(len(z_grid))
    for file_name in list_of_files:
        if 'stars' in file_name:
            data = np.genfromtxt(os.path.join(data_dir, file_name), dtype=dtype)
            for zz, cc in zip(data['z'], data['ct']):
                ix = int(np.round(zz/d_z))
                ct_grid[ix] += cc

    plt.figsize=(30,30)
    plt.plot(np.log10(z_grid), ct_grid)
    plt.xlabel('z (pc)')
    plt.ylabel('number of stars')
    plt.savefig(os.path.join(fig_dir, 'z_distribution_of_stars.png'))
