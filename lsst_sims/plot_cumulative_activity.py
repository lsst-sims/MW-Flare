"""
This script will take the fits fit_by_type.py and multiply by
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


if __name__ == "__main__":

    ct_dir = "query_results"
    gp_dir = "z_bin_fits"
    fig_dir = "plots"

    z_step = 25

    table_list = ['0870', '1100', '1160',
                  '1180','1200', '1220',
                  '1250', '1400']

    z_bin = []
    n_active_m4m6 = []
    n_total_m4m6 = []
    n_total = []
    type_ct = {}
    type_total = {}

    aa_dict = {}
    bb_dict = {}
    tau_dict = {}
    for spec_type in range(9):
        spec_name = 'M%d' % spec_type
        with open(os.path.join('type_fits', '%s_fit.txt' % spec_name)) as input_file:
            lines = input_file.readlines()
            params = lines[1].split()
            aa = float(params[0])
            tau = float(params[1])
            bb = float(params[2])

        aa_dict[spec_type] = aa
        bb_dict[spec_type] = bb
        tau_dict[spec_type] = tau

    for z_min in range(25, 201, z_step):

        z_bin.append(float(z_min)+0.5*float(z_step))
        n_active_m4m6.append(0.0)
        n_total_m4m6.append(0.0)
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
                        spec_class = 12
                    ct = int(vv[1])
                    if spec_class<9:
                        aa = aa_dict[spec_class]
                        bb = bb_dict[spec_class]
                        tau = tau_dict[spec_class]
                    else:
                        aa = aa_dict[8]
                        bb = bb_dict[8]
                        tau = tau_dict[8]
                    frac = aa*np.exp(-1.0*z_bin[-1]/tau) + bb

                    if spec_class>=4 and spec_class<=6:
                        n_active_m4m6[-1] += frac*ct
                        n_total_m4m6[-1] += ct
                    n_total[-1] += ct

                    if spec_class in type_ct:
                        type_ct[spec_class] += frac*ct
                        type_total[spec_class] += ct
                    else:
                        type_ct[spec_class] = frac*ct
                        type_total[spec_class] = ct

    z_bin = np.array(z_bin)
    n_active_m4m6 = np.array(n_active_m4m6)
    n_total_m4m6 = np.array(n_total_m4m6)
    n_total = np.array(n_total)
    total_active_m4m6 = n_active_m4m6.sum()
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
    hh, = plt.plot(z_bin, 0.9*np.cumsum(n_active_m4m6)/total_active_m4m6, marker='o', color='b')
    header_list.append(hh)
    title_list.append('this model')
    plt.xlabel('distance from Galactic plane (pc)', fontsize=10)
    plt.ylabel('cumulative active fraction', fontsize=10)
    plt.ylim(0,1.4)
    plt.xlim(0, 250)
    xticks = [xx for xx in range(0, 250, 10)]
    xlabels = ['%d' % xx if ii%10==0 else '' for ii, xx in enumerate(xticks)]
    plt.xticks(xticks, xlabels)
    yticks = [xx for xx in np.arange(0.0, 1.4, 0.1)]
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

    plt.subplot(2,2,3)
    data_dtype = np.dtype([('class', int), ('frac', float)])
    data = np.genfromtxt('data/active_fraction_by_type_West_et_al_2008.txt',
                         dtype=data_dtype)

    hh, = plt.plot(data['class'], data['frac'], color='r')
    header_list = [hh]
    label_list = ['West et al 2008 fig 3']

    type_arr = list(type_ct.keys())
    type_arr.sort()
    type_ct_arr = np.array([type_ct[cc] for cc in type_arr])
    type_total_arr = np.array([type_total[cc] for cc in type_arr])
    hh, = plt.plot(type_arr, type_ct_arr.astype(float)/type_total_arr, marker='o', linestyle='')
    header_list.append(hh)
    label_list.append('this model')
    hh, = plt.plot(type_arr, type_total_arr.astype(float)/type_total_arr.sum(), linestyle='--', color='k')
    header_list.append(hh)
    label_list.append('fraction of stars in type')
    plt.xlabel('spectral class', fontsize=10)
    plt.ylabel('fraction active', fontsize=10)
    xticks = range(10)
    xlabels = ['M%d' %ii for ii in xticks]
    xticks.append(12)
    xlabels.append('later')
    plt.xticks(xticks, xlabels,fontsize=10)
    plt.legend(header_list, label_list, fontsize=7, loc=2)
    plt.xlim(-2,13)

    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "hilton_2010_fig12.png"))
