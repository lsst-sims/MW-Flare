'''
A toy model for generating flares along a given line-of-sight in LSST

'''

import numpy as np
import pandas as pd


'''
Things I'll need:
- fields of view to compute over
    - use Trilegal fields for now
- Kepler flare light curve model (Davenport 2014)
- Kepler -> ugriz flare model or approximation
- LSST cadence, or cadence approximation
- Kepler-based SpT vs Flare Rate model (Davenport in prep)


'''

