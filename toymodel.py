'''
A toy model for generating flares along a given line-of-sight in LSST

'''


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


'''
Things I'll need:
- fields of view to compute over
    - use Trilegal fields for now
- Kepler flare light curve model (Davenport 2014)
- Kepler -> ugriz flare model or approximation
- LSST cadence, or cadence approximation
- Kepler-based SpT vs Flare Rate model (Davenport in prep)

'''



def generate_visits(Nvisits=900, tspan=10, stat=False):
    '''
    Use some very crude approximations for how visits will be spaced out:

    - Survey starts at midnight, time = 0.0
    - Can only observe at night, time > 0.75 | time < 0.25
    - Exposures are clustered around a season w/ a gaussian shape each year
    - Field is observable for first half of year, 0 < date < 182
    - On average, each field is hit every 3 days during observable season

    '''

    # generate random times for visit, between [0.75 and 0.25]
    time_of_day = np.random.random(Nvisits)/2. - 0.25

    date_of_year = np.floor(np.random.normal(loc=365./4., scale=365./7., size=Nvisits))

    year_of_obs = np.floor(np.random.random(Nvisits) * tspan) * 365.

    date_obs = time_of_day + date_of_year + year_of_obs

    date_obs.sort()

    if stat is True:
        print('mean time between visits:')
        print(np.mean(date_obs[1:] - date_obs[:-1]))

        print('median time between visits:')
        print(np.median(date_obs[1:] - date_obs[:-1]))

        plt.figure()
        _ = plt.hist(date_obs, bins=np.arange(date_obs.min(), date_obs.max(),7),
                     histtype='stepfilled', color='k')
        plt.xlabel('Time (days)')
        plt.ylabel('# Visits per Week')
        plt.show()

    return date_obs



