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
parser.add_argument('--duration', type=float, default=10.0,
                    help='the length (in years) of the light cuves '
                         'to simulate')
args = parser.parse_args()

if args.out_file is None:
    raise RuntimeError("You must specify out_file")

if not args.out_file.endswith('.npz'):
    args.out_file += '.npz'

t_start = time.time()

rng = np.random.RandomState(991)

cache = {}

for star_class in ('early_active', 'early_inactive',
                   'mid_active', 'mid_inactive', 'late_active'):

    for ix in range(4):
        tag = '%s_%d' % (star_class, ix)

        (cache['%s_time' % tag],
         cache['%s_u' % tag],
         cache['%s_g' % tag],
         cache['%s_r' % tag],
         cache['%s_i' % tag],
         cache['%s_z' % tag],
         cache['%s_y' % tag]) = light_curve_from_class(star_class, args.duration, rng)

        print(star_class,ix,(time.time()-t_start)/60.0)

with open(args.out_file, 'wb') as file_handle:
    np.savez(file_handle, **cache)
