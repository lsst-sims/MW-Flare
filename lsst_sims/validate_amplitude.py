"""
This script will try to recreate the amplitude-versus-energy and
amplitude-versus-duration panels from Hawley et al 2014 (ApJ 797, 12)
Figure 10
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mdwarf_utils import duration_from_energy, amplitude_from_duration_energy

import numpy as np

rng = np.random.RandomState(77)

n_flares = 10000
log_e = rng.random_sample(n_flares)*4.0 + 29.0
ee = np.power(10.0, log_e)

duration = duration_from_energy(ee, rng)
amplitude_u = amplitude_from_duration_energy(duration, ee)

log_ekp_quiescent = 30.67 # GJ1243 in Table 2 of Hawley et al 2014
ekp_quiescent = np.power(10.0, log_ekp_quiescent)

# Convert to Kepler amplitude.
# See paragraph before section 3 of Hawley et al 2014
amplitude_kp = amplitude_u/0.65 

plt.figsize = (30, 30)
plt.subplot(1,2,1)
plt.scatter(np.log10(ee/0.65), amplitude_kp/ekp_quiescent)
plt.xlabel('Log(E_kp) in ergs')
plt.ylabel('relative amplitude')
plt.yscale('log')
plt.xlim(29, 33)
plt.xticks(range(29,34))
plt.ylim(0.0005, 0.5)


plt.subplot(1,2,2)
plt.scatter(duration, amplitude_kp/ekp_quiescent)
plt.xlabel('duration (minutes)')
plt.ylabel('relative amplitude')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1,200)
plt.ylim(0.0005, 0.5)

plt.tight_layout()
plt.savefig('plots/amplitude_plot.png')
