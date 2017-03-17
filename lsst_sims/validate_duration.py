"""
This script will try to recreate the energy-versus-duration plot from
Hawley et al 2014 (ApJ 797, 12) Figure 10
"""

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from mdwarf_utils import duration_from_energy

import numpy as np

rng = np.random.RandomState(44)

n_flares = 10000
log_e = rng.random_sample(n_flares)*4.0 + 29.0
ee = np.power(10.0, log_e)

duration = duration_from_energy(ee, rng)

plt.figsize = (30, 30)

log_ekp = log_e - np.log10(0.65)

plt.scatter(log_ekp, duration, marker='o', color='k')

plt.yscale('log')

plt.xlabel('log(E_kepler) (ergs)')
plt.ylabel('duration (min)')
plt.xlim(29,33)
plt.xticks(range(29,33))
plt.ylim(1,200)

plt.savefig('plots/duration_plot.png')
