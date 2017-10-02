from __future__ import with_statement

import numpy as np
import os

from mdwarf_utils import light_curve_from_class
import time

t_start = time.time()

rng = np.random.RandomState(991)

cache = {}

for star_class in ('early_active', 'early_inactive',
                   'mid_active', 'mid_inactive', 'late_active'):

    for ix in range(4):
        tt, uu, gg, rr, ii, zz, yy = light_curve_from_class(star_class, 0.25, rng)

        tag = '%s_%d' % (star_class, ix)
        cache['%s_time' % tag] = tt
        cache['%s_u' % tag] = uu
        cache['%s_g' % tag] = gg
        cache['%s_r' % tag] = rr
        cache['%s_i' % tag] = ii
        cache['%s_z' % tag] = zz
        cache['%s_y' % tag] = yy

        print(star_class,ix,(time.time()-t_start)/60.0)


out_dir = os.path.join('/Users', 'danielsf', 'physics', 'lsst_160212',
                       'Development', 'sims_catUtils', 'workspace', 'mlt', 'data')
with open(os.path.join(out_dir, 'mdwarf_flare_light_curves_171002.npz'), 'wb') as file_handle:
    np.savez(file_handle, **cache)
