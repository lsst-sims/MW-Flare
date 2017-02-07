"""
This script will take the GP fits from fit_zbins.py and multiply by
the number of stars in each spectral class in each z bin to try to
reproduce the dashed line ("Active Stars") in Figure 12 of Hilton et al 2010
(AJ 140, 1402)
"""
from __future__ import with_statement

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os

from gaussian_process import ExpSquaredKernel, Covariogram
from gaussian_process import ForcedMeanGP

if __name__ == "__main__":

    ct_dir = "output"
    gp_dir = "z_bin_fits"
    fig_dir = "plots"

    z_step = 25

    table_list = ['0870', '1100', '1160',
                  '1180','1200', '1220',
                  '1250', '1400']

    z_bin = []
    n_active = []
    n_total = []

    kernel = ExpSquaredKernel(dim=1)
    covariogram = Covariogram(kernel)

    gp_dtype = np.dtype([('class', float), ('frac', float),
                         ('minus', float), ('plus', float)])

    for z_min in range(25, 201, z_step):
        gp_data_name = os.path.join("z_bins",
                                    "bin_%d_%d.txt" % (z_min, z_min+z_step))

        gp_data = np.genfromtxt(gp_data_name, dtype=gp_dtype)

        gp_name = os.path.join(gp_dir, "gp_%d_%d.txt" % (z_min, z_min+z_step))
        with open(gp_name, "r") as input_file:
            gp_lines = input_file.readlines()
            hp_line = gp_lines[1].split()[3:]
            hyper_params = np.array([float(hp_line[0]), float(hp_line[1])])
            nugget = []
            for line in gp_lines:
                if line[0] != '#':
                    break
                if 'nugget' in line:
                    nugget.append(float(line.split()[3]))

            nugget = np.array(nugget)
            for line in gp_lines:
                if line[0]!='#':
                    break
                if 'mean' in line:
                    forced_mean = float(line.split()[2])

            gp = ForcedMeanGP(covariogram, forced_mean)
            gp.covariogram.nugget = nugget
            gp.build(gp_data['class'], gp_data['frac'])
            z_bin.append(float(z_min)+0.5*float(z_step))
            n_active.append(0.0)
            n_total.append(0.0)

            for table in table_list:
                ct_name = os.path.join(ct_dir,
                                       'mdwarf_count_%s_%d_%d.txt' %
                                       (table, z_min, z_min+z_step))
                with open(ct_name, "r") as input_file:
                    for line in input_file:
                        vv = line.split()
                        if vv[0].startswith('M'):
                            spec_class = int(vv[0].replace('M','').replace(':',''))
                        else:
                            spec_class = 12.0
                        ct = int(vv[1])
                        frac = gp.regress([float(spec_class)])
                        n_active[-1] += frac[0]*ct
                        n_total[-1] += ct

    z_bin = np.array(z_bin)
    n_active = np.array(n_active)
    n_total = np.array(n_total)
    total_active = n_active.sum()
    plt.figsize = (30,30)

    control_dtype = np.dtype([('z', float),  ('frac', float)])
    control_data = np.genfromtxt('data/activity_rate_Hilton_et_al_2010.txt',
                                 dtype=control_dtype)

    plt.subplot(2,2,1)
    hh, = plt.plot(control_data['z'], control_data['frac'], color='r')
    header_list = [hh]
    title_list = ['Hilton et al. 2010']

    # mutliply by 0.9 because 0.9 of the active stars in
    # Hilton et al. Figure 12 occcur by the 225 pc mark,
    # where our data runs out
    hh, = plt.plot(z_bin, 0.9*np.cumsum(n_active)/total_active, marker='o', color='b')
    header_list.append(hh)
    title_list.append('this model')
    plt.xlabel('distance from Galactic plane (pc)', fontsize=10)
    plt.ylabel('cumulative active fraction', fontsize=10)
    plt.ylim(0,1.2)
    plt.xlim(0, 250)
    xticks = [xx for xx in range(0, 250, 10)]
    xlabels = ['%d' % xx if ii%10==0 else '' for ii, xx in enumerate(xticks)]
    plt.xticks(xticks, xlabels)
    yticks = [xx for xx in np.arange(0.0, 1.2, 0.1)]
    ylabels = ['%.1f' % xx if ii%4 ==0 else '' for ii, xx in enumerate(yticks)]
    plt.yticks(yticks, ylabels)
    plt.legend(header_list, title_list, fontsize=10)

    plt.subplot(2,2,2)
    plt.plot(z_bin, n_total, marker='x', color='k')
    plt.xlabel('distance from Galactic plane (pc)', fontsize=10)
    plt.ylabel('number of stars', fontsize=10)
    plt.xlim(0, 250)
    xticks = [xx for xx in range(0, 250, 10)]
    xlabels = ['%d' % xx if ii%10==0 else '' for ii, xx in enumerate(xticks)]
    plt.xticks(xticks, xlabels)
    #yticks = [xx for xx in np.arange(0.0, 1.1, 0.05)]
    #ylabels = ['%.1f' % xx if ii%4 ==0 else '' for ii, xx in enumerate(yticks)]
    #plt.yticks(yticks, ylabels)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "hilton_2010_fig12.png"))
