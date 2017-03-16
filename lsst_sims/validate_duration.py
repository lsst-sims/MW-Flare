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

d_min = np.interp(log_ekp,
                  duration_from_energy._model_data['ekp'],
                  duration_from_energy._model_data['min'])

d_max = np.interp(log_ekp,
                  duration_from_energy._model_data['ekp'],
                  duration_from_energy._model_data['max'])

d_mean = np.interp(log_ekp,
                   duration_from_energy._model_data['ekp'],
                   duration_from_energy._model_data['mean'])

plt.plot(log_ekp, d_min)
plt.plot(log_ekp, d_mean, linestyle='--')
plt.plot(log_ekp, d_max)


plt.yscale('log')

plt.xlabel('log(E_kepler) (ergs)')
plt.ylabel('duration (min)')

plt.savefig('plots/duration_plot.png')
