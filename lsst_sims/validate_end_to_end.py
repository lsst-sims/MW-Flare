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

def make_distribution(xx_in, dd, color='k'):
    xx_dex = np.round(xx_in/dd).astype(int)
    unq, ct = np.unique(xx_dex, return_counts=True)

    sorted_dex = np.argsort(unq)

    return plt.plot(unq[sorted_dex]*dd, ct[sorted_dex], color=color)


def make_density_plot(xx_in, yy_in, dd):
    cmin = 0
    cmax = 150
    dc = 25

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

    xx_grid, yy_grid = np.mgrid[slice(xx_arr.min()*dd, xx_arr.max()*dd+dd, dd),
                                slice(yy_arr.min()*dd, yy_arr.max()*dd+dd, dd)]

    cc_grid = np.ones(xx_grid.shape)*(-99.0)

    for ix, iy, cc in zip(xx_arr, yy_arr, ct_arr):
        cc_grid[ix-xx_arr.min()][iy-yy_arr.min()] = cc

    cc_masked = np.ma.masked_values(cc_grid, -99.0)
    plt.pcolormesh(xx_grid, yy_grid, cc_masked,
                   cmap=plt.cm.gist_ncar,
                   edgecolor='')

    cb = plt.colorbar()
    plt.clim(cmin,cmax)
    cb.set_ticks(np.arange(cmin,cmax+10,dc))

    return None

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
simulation_length = control_data['stop'].max()-control_data['start'].min()

t_flare, e_flare = draw_energies('mid_active', simulation_length, rng)

duration = duration_from_energy(e_flare, rng)

amp = amplitude_from_duration_energy(duration, e_flare)

ekp_flare = e_flare/0.65 # paragraph before section 3 of Hawley et al 2014
amp_kp = amp/0.65

log_ekp_flare = np.log10(ekp_flare)

amp_rel = amp_kp/np.power(10.0,log_ekp_quiescent)

plt.figsize = (30,30)

dx = 0.05

plt.subplot(3,2,1)

make_density_plot(log_ekp_flare, np.log10(amp_rel), dx)

plt.title('simulation', fontsize=10)
plt.xlabel('Log(E_kp) in ergs',fontsize=10)
plt.ylabel('Log(relative amplitude)', fontsize=10)
plt.ylim(-4, 1)
plt.yticks(range(-4, 1))
plt.xlim(29,35)
plt.xticks(range(29,35))


plt.subplot(3,2,2)
make_density_plot(control_log_ekp, np.log10(control_amp), dx)

plt.title('data', fontsize=10)
plt.xlabel('Log(E_kp) in ergs',fontsize=10)
plt.ylabel('Log(relative amplitude)', fontsize=10)
plt.ylim(-4, 1)
plt.yticks(range(-4, 1))
plt.xlim(29,35)
plt.xticks(range(29,35))


plt.subplot(3,2,3)
make_density_plot(log_ekp_flare, np.log10(duration), dx)

plt.title('simulation', fontsize=10)
plt.xlabel('Log(E_kp) in ergs',fontsize=10)
plt.ylabel('Log(duration) (minutes)', fontsize=10)
plt.xlim(29,35)
plt.ylim(-1,3)
plt.yticks(range(-1,3))
plt.xticks(range(29,35))

plt.subplot(3,2,4)
make_density_plot(control_log_ekp, np.log10(control_duration), dx)

plt.title('data', fontsize=10)
plt.xlabel('Log(E_kp) in ergs', fontsize=10)
plt.ylabel('Log(duration) (minutes)', fontsize=10)
plt.xlim(29,35)
plt.ylim(-1,3)
plt.yticks(range(-1,3))
plt.xticks(range(29,35))

plt.subplot(3,2,5)
make_density_plot(np.log10(duration), np.log10(amp_rel), dx)

plt.title('simulation',fontsize=10)
plt.xlabel('Log(duration) (minutes)', fontsize=10)
plt.ylabel('Log(relative amplitude)', fontsize=10)
plt.ylim(-4, 0)
plt.yticks(range(-4,0))
plt.xlim(0,2.5)
plt.xticks(range(0,3))


plt.subplot(3,2,6)
make_density_plot(np.log10(control_duration), np.log10(control_amp),dx)

plt.title('data',fontsize=10)
plt.xlabel('Log(duration) (minutes)', fontsize=10)
plt.ylabel('Log(relative amplitude)', fontsize=10)
plt.ylim(-4, 0)
plt.yticks(range(-4,0))
plt.xlim(0,2.5)
plt.xticks(range(0,3))


plt.tight_layout()
plt.savefig('plots/end_to_end_2D.png')
print 'control ',len(control_data),' sim ',len(amp_rel)
plt.close()

plt.figsize = (30, 30)
plt.subplot(2,2,1)
header_list = []
label_list = []
hh, = make_distribution(log_ekp_flare, 0.1, color='k')
header_list.append(hh)
label_list.append('simulations')
hh, = make_distribution(control_log_ekp, 0.1, color='r')
header_list.append(hh)
label_list.append('data')
plt.xlabel('Log(E_kp)')

plt.legend(header_list, label_list, loc=0)

plt.subplot(2,2,2)
make_distribution(np.log10(duration), 0.1, color='k')
make_distribution(np.log10(control_duration), 0.1, color='r')
plt.xlabel('Log(duration)')

plt.subplot(2,2,3)
make_distribution(np.log10(amp_rel), 0.1, color='k')
make_distribution(np.log10(control_amp), 0.1, color='r')
plt.xlabel('Log(amplitude)')

plt.tight_layout()
plt.savefig('plots/end_to_end_1D.png')
plt.close()
