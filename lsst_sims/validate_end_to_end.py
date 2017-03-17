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

def make_density_data(xx_in, yy_in, dd):
    data_x = np.round(xx_in/dd).astype(int)
    data_y = np.round(yy_in/dd).astype(int)

    xmax = data_x.max()
    ymax = data_y.max()
    xmin = data_x.min()
    ymin = data_y.min()
    factor = int(np.power(10.0, np.round(np.log10(ymax-ymin)+1.0)))
    dex_arr = (data_x-xmin)*factor + data_y-ymin

    unq, counts = np.unique(dex_arr, return_counts=True)

    x_unq = xmin + unq//factor
    y_unq = ymin + unq % factor

    grid = {}

    for xx, yy, cc in zip(x_unq, y_unq, counts):
        if xx not in grid:
            grid[xx]= {}

        if yy not in grid[xx]:
            grid[xx][yy] = cc
        else:
            grid[xx][yy] += cc

    xx_arr = []
    yy_arr = []
    ct_arr = []
    for xx in grid:
        for yy in grid[xx]:
            xx_arr.append(xx)
            yy_arr.append(yy)
            ct_arr.append(grid[xx][yy])

    ct_arr = np.array(ct_arr)
    xx_arr = np.array(xx_arr)
    yy_arr = np.array(yy_arr)

    sorted_dex = np.argsort(ct_arr)

    xx_arr = xx_arr[sorted_dex]*dd
    yy_arr = yy_arr[sorted_dex]*dd
    ct_arr = ct_arr[sorted_dex]
    ct_arr = np.cumsum(ct_arr)/float(ct_arr.sum())

    return xx_arr, yy_arr, ct_arr


cmin = 0
cmax = 1
dc = 0.2


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

control_duration = (control_data['stop']-control_data['start'])*24.0*60.0
control_log_ekp = log_ekp_quiescent + np.log10(control_data['e_dur'])
control_amp = control_data['amp']

rng = np.random.RandomState(813)
t_flare, e_flare = draw_energies('mid_active', 720.0, rng)

duration = duration_from_energy(e_flare, rng)

amp = amplitude_from_duration_energy(duration, e_flare)

ekp_flare = e_flare/0.65 # paragraph before section 3 of Hawley et al 2014
amp_kp = amp/0.65

log_ekp_flare = np.log10(ekp_flare)

amp_rel = amp_kp/np.power(10.0,log_ekp_quiescent)

plt.figsize = (30,30)

dx = 0.05

plt.subplot(3,2,1)

xx, yy, cc = make_density_data(log_ekp_flare, np.log10(amp_rel), dx)

plt.scatter(xx, yy, c=cc, edgecolor='',
            cmap=plt.cm.gist_ncar,
            s=5)
cb = plt.colorbar()
plt.clim(cmin,cmax)
cb.set_ticks(np.arange(cmin,cmax+10,dc))
plt.title('simulation', fontsize=10)
plt.xlabel('Log(E_kp) in ergs',fontsize=10)
plt.ylabel('Log(relative amplitude)', fontsize=10)
plt.ylim(-4, 1)
plt.yticks(range(-4, 1))
plt.xlim(29,35)
plt.xticks(range(29,35))


plt.subplot(3,2,2)
xx, yy, cc = make_density_data(control_log_ekp, np.log10(control_amp), dx)

plt.scatter(xx, yy, c=cc, edgecolor='',
            cmap=plt.cm.gist_ncar,
            s=5)
cb = plt.colorbar()
plt.clim(cmin,cmax)
cb.set_ticks(np.arange(cmin,cmax+10,dc))
plt.title('data', fontsize=10)
plt.xlabel('Log(E_kp) in ergs',fontsize=10)
plt.ylabel('Log(relative amplitude)', fontsize=10)
plt.ylim(-4, 1)
plt.yticks(range(-4, 1))
plt.xlim(29,35)
plt.xticks(range(29,35))


plt.subplot(3,2,3)
xx, yy, cc = make_density_data(log_ekp_flare, np.log10(duration), dx)

plt.scatter(xx, yy, c=cc, edgecolor='',
            cmap=plt.cm.gist_ncar,
            s=5)
cb = plt.colorbar()
plt.clim(cmin,cmax)
cb.set_ticks(np.arange(cmin,cmax+10,dc))
plt.title('simulation', fontsize=10)
plt.xlabel('Log(E_kp) in ergs',fontsize=10)
plt.ylabel('Log(duration) (minutes)', fontsize=10)
plt.xlim(29,35)
plt.ylim(-1,3)
plt.yticks(range(-1,3))
plt.xticks(range(29,35))

plt.subplot(3,2,4)

xx, yy, cc = make_density_data(control_log_ekp, np.log10(control_duration), dx)

plt.scatter(xx, yy, c=cc, edgecolor='',
            cmap=plt.cm.gist_ncar,
            s=5)
cb = plt.colorbar()
plt.clim(cmin,cmax)
cb.set_ticks(np.arange(cmin,cmax+10,dc))
plt.title('data', fontsize=10)
plt.xlabel('Log(E_kp) in ergs', fontsize=10)
plt.ylabel('Log(duration) (minutes)', fontsize=10)
plt.xlim(29,35)
plt.ylim(-1,3)
plt.yticks(range(-1,3))
plt.xticks(range(29,35))

plt.subplot(3,2,5)
xx, yy, cc = make_density_data(np.log10(duration), np.log10(amp_rel), dx)

plt.scatter(xx, yy, c=cc, edgecolor='',
            cmap=plt.cm.gist_ncar,
            s=10)
cb = plt.colorbar()
plt.clim(cmin,cmax)
cb.set_ticks(np.arange(cmin,cmax+10,dc))
plt.title('simulation',fontsize=10)
plt.xlabel('Log(duration) (minutes)', fontsize=10)
plt.ylabel('Log(relative amplitude)', fontsize=10)
plt.ylim(-4, 0)
plt.yticks(range(-4,0))
plt.xlim(0,2.5)
plt.xticks(range(0,3))


plt.subplot(3,2,6)
xx, yy, cc = make_density_data(np.log10(control_duration), np.log10(control_amp),
                               dx)

plt.scatter(xx, yy, c=cc, edgecolor='',
            cmap=plt.cm.gist_ncar,
            s=10)
cb = plt.colorbar()
plt.clim(cmin,cmax)
cb.set_ticks(np.arange(cmin,cmax+10,dc))
plt.title('data',fontsize=10)
plt.xlabel('Log(duration) (minutes)', fontsize=10)
plt.ylabel('Log(relative amplitude)', fontsize=10)
plt.ylim(-4, 0)
plt.yticks(range(-4,0))
plt.xlim(0,2.5)
plt.xticks(range(0,3))


plt.tight_layout()
plt.savefig('plots/end_to_end_plot.png')
print 'control ',len(control_data),' sim ',len(amp_rel)
