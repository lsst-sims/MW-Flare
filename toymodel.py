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

# flare morphology, taken from "appaloosa" (Davenport 2014)
import aflare

from LSSToy import generate_visits


def run_field(file):
    '''
    for this TRILEGAL field:
    - generate a cadence model
    - for every star generate simulated flares as a function of color (mass) and age, based on in prep work

    '''

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
    all_fields()
