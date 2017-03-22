"""
This script will try to recreate the energy-versus-duration plot from
Hawley et al 2014 (ApJ 797, 12) Figure 10
"""

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from mdwarf_utils import duration_from_energy
from plot_utils import make_density_plot, make_distribution_plot

import numpy as np

dtype = np.dtype([('start_dex', float), ('stop_dex', float),
                  ('start', float), ('stop', float), ('peak', float),
                  ('rise', float), ('decay', float), ('amp', float),
                  ('e_dur', float), ('e_dur_rise', float),
                  ('e_dur_decay', float), ('flag', float),
                  ('ppl_flare', float), ('ppl_month', float),
                  ('components', float)])

control_data = np.genfromtxt('data/gj1243_master_flares.tbl', dtype=dtype)
valid = np.where(np.logical_not(np.isnan(control_data['amp'])))
control_data = control_data[valid]
valid = np.where(np.logical_and(control_data['amp']>0.0, control_data['e_dur']>0.0))
control_data = control_data[valid]

log_ekp_quiescent = 30.67
control_log_ekp = log_ekp_quiescent + np.log10(control_data['e_dur'])

log_eu = control_log_ekp + np.log10(0.65)

eu = np.power(10.0, log_eu)
rng = np.random.RandomState(88)

duration = duration_from_energy(eu, rng)

dx = 0.05
plt.figsize = (30, 30)
plt.subplot(2,2,1)
make_density_plot(log_eu, np.log10(duration), dx)
plt.xlabel('Log(E_U)')
plt.ylabel('Log(simulated duration)')
plt.xlim(26,34)
plt.ylim(-0.5, 2.5)

plt.subplot(2,2,2)
control_duration = 24.0*60.0*(control_data['stop']-control_data['start'])
make_density_plot(log_eu, np.log10(control_duration), dx)
plt.xlabel('Log(E_U)')
plt.ylabel('Log(duration)')
plt.xlim(26,34)
plt.ylim(-0.5, 2.5)


plt.subplot(2,2,3)
hh, = make_distribution_plot(np.log10(duration), dx, color='k')
header_list = [hh]
label_list = ['simulated']
hh, = make_distribution_plot(np.log10(control_duration), dx, color='r')
header_list.append(hh)
label_list.append('data')
plt.legend(header_list, label_list, loc=0)
plt.xlabel('Log(duration in minutes)')
plt.ylabel('N')

plt.subplot(2,2,4)
make_distribution_plot(np.log10(np.abs(duration-control_duration)), dx)
plt.xlabel('Log( |duration - simulated duration| )')
plt.ylabel('N')

plt.tight_layout()
plt.savefig('plots/duration_plot.png')
