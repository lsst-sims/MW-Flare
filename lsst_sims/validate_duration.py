import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

from mdwarf_utils import duration_from_energy

import numpy as np

rng = np.random.RandomState()

n_flares = 10000
log_e = rng.random_sample(n_flares)*4.0 + 29.0
ee = np.power(10.0, log_e)

duration = duration_from_energy(ee, rng)

plt.figsize = (30, 30)

log_ekp = log_e - np.log10(0.65)

plt.scatter(log_ekp, duration, marker='o', color='k')

log_e = np.sort(log_e)
log_ekp = np.sort(log_ekp)

mm = duration_from_energy._models['min']['m']
bb = duration_from_energy._models['min']['b']

d_min = np.power(10.0, mm*log_e + bb)

mm = duration_from_energy._models['max']['m']
bb = duration_from_energy._models['max']['b']

d_max = np.power(10.0, mm*log_e + bb)

mm = duration_from_energy._models['mean']['m']
bb = duration_from_energy._models['mean']['b']

d_mean = np.power(10.0, mm*log_e + bb)

plt.plot(log_ekp, d_min)
plt.plot(log_ekp, d_mean, linestyle='--')
plt.plot(log_ekp, d_max)


plt.yscale('log')

plt.xlabel('log(E_kepler) (ergs)')
plt.ylabel('duration (min)')

plt.savefig('plots/duration_plot.png')
