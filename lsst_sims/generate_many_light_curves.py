from __future__ import with_statement

import numpy as np
import os

from mdwarf_utils import light_curve_from_class
import time

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--out_file', type=str, default=None,
                    help='the name of the file in which to store '
                         'the output')
args = parser.parse_args()

if args.out_file is None:
    raise RuntimeError("You must specify out_file")

if not args.out_file.endswith('.npz'):
    args.out_file += '.npz'

t_start = time.time()

rng = np.random.RandomState(991)

cache = {}

duration = 10.0 # in years

for star_class in ('early_active', 'early_inactive',
                   'mid_active', 'mid_inactive', 'late_active'):

    for ix in range(4):
        tt, uu, gg, rr, ii, zz, yy = light_curve_from_class(star_class, duration, rng)

        tag = '%s_%d' % (star_class, ix)
        cache['%s_time' % tag] = tt
        cache['%s_u' % tag] = uu
        cache['%s_g' % tag] = gg
        cache['%s_r' % tag] = rr
        cache['%s_i' % tag] = ii
        cache['%s_z' % tag] = zz
        cache['%s_y' % tag] = yy

        print(star_class,ix,(time.time()-t_start)/60.0)

with open(args.out_file, 'wb') as file_handle:
    np.savez(file_handle, **cache)
