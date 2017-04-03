"""
This script will draw the energies of flares for different classes of flaring
stars and plot the resulting cumulative distribution against the fits from
Table 4.3 of Eric Hilton's PhD disseration.  The plot will be

plots/energy_dist_validation.png
"""

from __future__ import with_statement

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import os
import numpy as np

from mdwarf_utils import draw_energies

class_list = []
alpha_list = []
beta_list = []
logemin_list = []
logemax_list = []
with open(os.path.join('data', 'hilton_phd_table_4.3.txt'), 'r') as param_file:
    input_lines = param_file.readlines()
    for line in input_lines:
        if line[0] != '#':
            vv = line.split()
            class_list.append(vv[0])
            alpha_list.append(float(vv[1]))
            beta_list.append(float(vv[2]))
            logemin_list.append(float(vv[3]))
            logemax_list.append(float(vv[4]))

rng = np.random.RandomState(118)

duration = 3650.0

color_list = ['k', 'r', 'b', 'g', 'c']

plt.figsize = (30, 30)
header_list = []
label_list = []
for ic, color in zip(range(len(class_list)), color_list):

    star_class = class_list[ic]
    alpha = alpha_list[ic]
    beta = beta_list[ic]
    logemin = logemin_list[ic]
    logemax = logemax_list[ic]

    tt, ee = draw_energies(star_class, duration, rng)
    print '\n',star_class
    print 'log(E_min/max): ',np.log10(ee.min()),np.log10(ee.max())
    print 'should be: ',logemin, logemax, logemin+0.1*(logemax-logemin)

    log_ee_control = np.arange(logemin, logemax, 0.1*(logemax-logemin))
    cum_ct = []
    for log_ee in log_ee_control:
        valid = np.where(ee>np.power(10.0, log_ee))
        n_valid = len(valid[0])
        cum_ct.append(float(n_valid)/(duration*24.0))

    cum_ct = np.array(cum_ct)

    hh, = plt.plot(log_ee_control, alpha+beta*log_ee_control, color=color)
    header_list.append(hh)
    label_list.append(star_class.replace('_',' '))
    plt.scatter(log_ee_control, np.log10(cum_ct), color=color, marker='+',
                s=55)

    overlarge = np.where(ee>np.power(10.0, 33.9))
    print 'number of flares set to e_abs_max: ',len(overlarge[0]), \
           ' out of ',len(ee)

    above_max = np.where(ee>np.power(10.0, logemax))
    print 'number of flares greater than emax: ',len(above_max[0]), \
          ' out of ',len(ee)

plt.legend(header_list, label_list, loc=0)
plt.xlabel('Log(Johnson U Band energy in ergs')
plt.ylabel('Log(cumulative number of flares)')
plt.text(27.1, -2.7, 'Markers are simulation results.\n'
         'Lines are the distributions being simulated.')
plt.savefig(os.path.join('plots', 'energy_dist_validation.png'))
