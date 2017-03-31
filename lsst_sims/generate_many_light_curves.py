from __future__ import with_statement

import numpy as np

from mdwarf_utils import light_curve_from_class
import time

t_start = time.time()

rng = np.random.RandomState(991)

tt_0 = None

cache = {}

for star_class in ('early_active', 'early_inactive',
                   'mid_active', 'mid_inactive', 'late_active'):

    for ix in range(3):
        tt, uu, gg, rr, ii, zz, yy, t_peak = light_curve_from_class(star_class, 5.0, rng)

        tag = '%s_%d' % (star_class, ix)
        cache['%s_time' % tag] = tt
        cache['%s_u' % tag] = uu
        cache['%s_g' % tag] = gg
        cache['%s_r' % tag] = rr
        cache['%s_i' % tag] = ii
        cache['%s_z' % tag] = zz
        cache['%s_y' % tag] = yy

        print(star_class,ix,(time.time()-t_start)/60.0)


with open('output/npz_light_curves.npz', 'wb') as file_handle:
    np.savez(file_handle, **cache)
