from __future__ import with_statement

import numpy as np

from mdwarf_utils import light_curve_from_class
import time

t_start = time.time()

rng = np.random.RandomState(991)

cache = {}

for star_class in ('early_active', 'early_inactive',
                   'mid_active', 'mid_inactive', 'late_active'):

    for ix in range(4):
        tt, uu, gg, rr, ii, zz, yy, t_peak = light_curve_from_class(star_class, 10.0, rng)

        tag = '%s_%d' % (star_class, ix)
        cache['%s_time' % tag] = tt
        cache['%s_u' % tag] = uu
        cache['%s_g' % tag] = gg
        cache['%s_r' % tag] = rr
        cache['%s_i' % tag] = ii
        cache['%s_z' % tag] = zz
        cache['%s_y' % tag] = yy

        print(star_class,ix,(time.time()-t_start)/60.0)


with open('/astro/store/pogo4/danielsf/mlt_flares/mdwarf_flare_light_curves_170927.npz', 'wb') as file_handle:
    np.savez(file_handle, **cache)
