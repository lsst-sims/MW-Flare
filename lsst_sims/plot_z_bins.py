import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os

if __name__ == "__main__":

    plot_dir = "plots"
    fit_dir = "z_bin_fits"
    data_dir = "z_bins"

    data_dtype = np.dtype([('class', float), ('frac', float),
                           ('lower', float), ('upper', float)])

    fit_dtype = np.dtype([('class', float), ('frac', float)])

    plt.figsize = (30,30)
    for i_fig, z_min in enumerate(range(0, 201, 25)):
        plt.subplot(3,3,i_fig+1)
        data_name = os.path.join(data_dir, "bin_%d_%d.txt" % (z_min, z_min+25))
        data = np.genfromtxt(data_name, dtype=data_dtype)

        fit_name = os.path.join(fit_dir, "gp_%d_%d.txt" % (z_min, z_min+25))
        fit = np.genfromtxt(fit_name, dtype=fit_dtype)

        plt.errorbar(data['class'], data['frac'],
                     yerr=np.array([data['frac']-data['lower'],
                                    data['upper']-data['frac']]),
                     linestyle='', marker='o', color='k')

        plt.plot(fit['class'], fit['frac'], color='r')

        plt.title('%d < |z| < %d pc' % (z_min, z_min+25), fontsize=10)
        plt.xlim(-1,13)
        plt.ylim(0,1.1)
        if i_fig==0:
            plt.xlabel('M dwarf class', fontsize=10)
            plt.ylabel('fraction spectrally active', fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'west_et_al_bins.png'))
