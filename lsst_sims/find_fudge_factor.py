"""
According to an email from Jim Davenport:

------
If I remember correctly what I said, it's that spectroscopic
activity doesn't == flaring activity. So stars could be "magnetically active"
as determined by Halpha emission in their spectra, but not "flare active"
(not as well defined)

The evidence of this is Hilton et al (2010),
http://adsabs.harvard.edu/abs/2010AJ....140.1402H (Fig 12) where they
find the scale height of flare stars is lower than for "active" stars.

Of course as the specialist I have to say "this is probably more
complicated"...

So my off-the-cuff advice:
do exactly as you propose, use West et al. to make "active vs inactive"
as a function of height, assign flares, etc. BUT, adopt some correction
(or ignorance) scaling factor that makes the scale height lower. I don't
know what that # should be, but you could extract something plausible
out of Fig 12 from Hilton2010 (e.g. 80pc instead of 120pc)
------

This script will take the fits to West et al. generated by
fit_by_type.py and try to determine the factor by which to multiply each
spectral type's scale height (tau) to recreate the "flare stars" curve in
Figure 12 of Hilton et al. 2010.

"""

from __future__ import with_statement

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
import numpy as np

from fit_activity_level import _find_fudge_factor


if __name__ == "__main__":

    data_dir = "data"
    plot_dir = "plots"

    dtype = np.dtype([('z', float), ('frac', float)])

    active_data = np.genfromtxt(os.path.join(data_dir,
                                             'activity_rate_Hilton_et_al_2010.txt'),
                                dtype=dtype)

    flare_data = np.genfromtxt(os.path.join(data_dir,
                                            'flare_rate_Hilton_et_al_2010.txt'),
                               dtype=dtype)

    (fudge_factor,
     tau_flare, offset_flare,
     tau_active, offset_active) = _find_fudge_factor()

    with open(os.path.join('type_fits','fudge_factor.txt'), 'w') as out_file:
        out_file.write('# factor to multiply tau by to get "flaring active" stars\n')
        out_file.write('%.9g\n' % (tau_flare/tau_active))

    plt.figsize = (30,30)
    hh, = plt.plot(active_data['z'], active_data['frac'], color='k')
    header_list = [hh]
    label_list = ['active stars']
    plt.plot(active_data['z'], 1.0-np.exp(-1.0*(active_data['z']-offset_active)/tau_active),
             color='k', linestyle='--')

    hh, = plt.plot(flare_data['z'], flare_data['frac'], color='r')
    header_list.append(hh)
    label_list.append('flaring stars')
    plt.plot(flare_data['z'], 1.0-np.exp(-1.0*(flare_data['z']-offset_flare)/tau_flare),
             color='r', linestyle='--')

    plt.xlabel('z(pc)')
    plt.ylabel('cumulative fraction')
    plt.legend(header_list, label_list, loc=0)
    plt.savefig(os.path.join('plots', 'fit_fudge_factor.png'))
