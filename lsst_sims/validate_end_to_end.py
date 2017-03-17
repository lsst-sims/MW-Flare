"""
This script will try to recreate all 3 panels of Figure 10 in
Hawley et al 2014 (ApJ 797, 121) by simulating one 'mid_active' star
for 360 days and plotting the characteristics of the resulting
flares
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

from mdwarf_utils import (draw_energies, duration_from_energy,
                          amplitude_from_duration_energy)

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

log_ekp_quiescent = 30.67

control_duration = (control_data['stop']-control_data['start'])*24.0*60.0
control_log_ekp = log_ekp_quiescent + np.log10(control_data['e_dur'])
control_amp = control_data['amp']

rng = np.random.RandomState(813)
t_flare, e_flare = draw_energies('mid_active', 360.0, rng)

duration = duration_from_energy(e_flare, rng)

amp = amplitude_from_duration_energy(duration, e_flare)

ekp_flare = e_flare/0.65 # paragraph before section 3 of Hawley et al 2014
amp_kp = amp/0.65

log_ekp_flare = np.log10(ekp_flare)

amp_rel = amp_kp/np.power(10.0,log_ekp_quiescent)

plt.figsize = (30,30)

plt.subplot(3,2,1)
plt.scatter(log_ekp_flare, amp_rel, color='b')
plt.xlabel('Log(E_kp) in ergs')
plt.ylabel('relative amplitude')
plt.yscale('log')
plt.ylim(0.0005, amp_rel.max())
plt.xlim(29,35)
plt.xticks(range(29,35))

plt.subplot(3,2,2)
plt.scatter(control_log_ekp, control_amp, color='r')
plt.xlabel('Log(E_kp) in ergs')
plt.ylabel('relative amplitude')
plt.yscale('log')
plt.ylim(0.0005, amp_rel.max())
plt.xlim(29,35)
plt.xticks(range(29,35))

plt.subplot(3,2,3)
plt.scatter(log_ekp_flare, duration, color='b')
plt.xlabel('Log(E_kp) in ergs')
plt.ylabel('duration (minutes)')
plt.yscale('log')
plt.xlim(29,35)
plt.ylim(0.1,1000.0)
plt.xticks(range(29,35))

plt.subplot(3,2,4)
plt.scatter(control_log_ekp, control_duration, color='r')
plt.xlabel('Log(E_kp) in ergs')
plt.ylabel('duration (minutes)')
plt.yscale('log')
plt.xlim(29,35)
plt.ylim(0.1,1000.0)
plt.xticks(range(29,35))

plt.subplot(3,2,5)
plt.scatter(duration, amp_rel, color='b')
plt.xlabel('duration (minutes)')
plt.ylabel('relative amplitude')
plt.yscale('log')
plt.xscale('log')
plt.ylim(0.0001, 0.2)
plt.xlim(1, 1000.0)

plt.subplot(3,2,6)
plt.scatter(control_duration, control_amp, color='r')
plt.xlabel('duration (minutes)')
plt.ylabel('relative amplitude')
plt.yscale('log')
plt.xscale('log')
plt.ylim(0.0001, 0.2)
plt.xlim(1, 1000.0)

plt.tight_layout()
plt.savefig('plots/end_to_end_plot.png')
print 'control ',len(control_data),' sim ',len(amp_rel)
