'''
A toy model for generating flares along a given line-of-sight in LSST

Things I'll need:
- fields of view to compute over
    - use Trilegal fields for now
- Kepler flare light curve model (Davenport 2014)
- Kepler -> ugriz flare model or approximation
- LSST cadence, or cadence approximation
- Kepler-based SpT vs Flare Rate model (Davenport in prep)

'''


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import LSSToy
import flare_prob


def downsample(time,flux):
    '''
    take super-sampled LC (from flare_prob), uses simple linear interpretation
    to down-sample to LSST cadence.

    Assumes 10 years of LSST with 900 visits
    '''

    tout = LSSToy.generate_visits()

    fout = np.interp(tout, time, flux)

    return tout, fout


def run_field(file, ):
    '''
    for this TRILEGAL field:
    - generate a cadence model
    - for every star generate simulated flares as a function of color (mass) and age, based on in prep work

    '''

    if file is 'test':
        print('Doing a sweep of alpha values [0.01, 0.1, 1]')

        alpha = [0.01, 0.1, 1.0]
        traw, fraw = flare_prob.SuperLC(dur=0.1, repeat=100, ffd_alpha=0.1)
        time, flux = downsample(traw, fraw)

    else:
        df = pd.read_table(file)

    return


def all_fields(models='index.txt'):
    dir = 'trilegal_models/'
    files = np.loadtxt(dir+models, comments='#', unpack=True, usecols=(0,), delimiter=',', dtype=np.str)

    for k in range(len(files)):
        print('running '+ dir + files[k])
        run_field(dir+files[k])

    return


if __name__ == "__main__":
    # all_fields()
    run_field('test')
