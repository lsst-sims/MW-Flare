from __future__ import with_statement

import numpy as np
import os

from mdwarf_utils import light_curve_from_params, generate_light_curve_params
from mdwarf_utils import light_curve_from_class
from mdwarf_utils import get_clean_dexes
import time

import multiprocessing as mproc

import argparse


def cache_lc(t_peak_arr=None, fwhm_arr=None, amp_arr=None,
             tag=None, out_dict=None, lock=None):

    (tt, uu,
     gg, rr,
     ii, zz,
     yy) = light_curve_from_params(t_peak_arr, fwhm_arr, amp_arr)

    out_dexes = get_clean_dexes(tt, gg, t_peak_arr, fwhm_arr)

    lock.acquire()
    out_dict['%s_time' % tag] = tt[out_dexes]
    out_dict['%s_u' % tag] = uu[out_dexes]
    out_dict['%s_g' % tag] = gg[out_dexes]
    out_dict['%s_r' % tag] = rr[out_dexes]
    out_dict['%s_i' % tag] = ii[out_dexes]
    out_dict['%s_z' % tag] = zz[out_dexes]
    out_dict['%s_y' % tag] = yy[out_dexes]
    lock.release()

parser = argparse.ArgumentParser()

parser.add_argument('--out_file', type=str, default=None,
                    help='the name of the file in which to store '
                         'the output')

parser.add_argument('--duration', type=float, default=10.0,
                    help='the length (in years) of the light cuves '
                         'to simulate')

parser.add_argument('--n_proc', type=int, default=1,
                    help='the number of independent processes to spawn')

args = parser.parse_args()

if args.out_file is None:
    raise RuntimeError("You must specify out_file")

if not args.out_file.endswith('.npz'):
    args.out_file += '.npz'

t_start = time.time()

rng = np.random.RandomState(991)

manager = mproc.Manager()

cache = manager.dict()

lock = mproc.Lock()
process_list = []

for star_class in ('late_active', 'early_inactive', 'early_active',
                   'mid_inactive', 'mid_active'):

    for ix in range(4):
        tag = '%s_%d' % (star_class, ix)

        if args.n_proc == 1:
            (t_peak_arr,
             fwhm_arr,
             amplitude_arr) = generate_light_curve_params(star_class, args.duration, rng)

            (tt,
             uu, gg, rr,
             ii, zz, yy) = light_curve_from_params(t_peak_arr, fwhm_arr, amplitude_arr)


            clean_dexes = get_clean_dexes(tt, gg, t_peak_arr, fwhm_arr)
            tt = tt[clean_dexes]
            uu = uu[clean_dexes]
            gg = gg[clean_dexes]
            rr = rr[clean_dexes]
            ii = ii[clean_dexes]
            zz = zz[clean_dexes]
            yy = yy[clean_dexes]

            cache['%s_time' % tag] = tt
            cache['%s_u' % tag] = uu
            cache['%s_g' % tag] = gg
            cache['%s_r' % tag] = rr
            cache['%s_r' % tag] = ii
            cache['%s_r' % tag] = zz
            cache['%s_z' % tag] = yy

        else:
            (t_peak_arr,
             fwhm_arr,
             amplitude_arr) = generate_light_curve_params(star_class, args.duration, rng)


            p = mproc.Process(target=cache_lc,
                              kwargs={'t_peak_arr': t_peak_arr,
                                      'fwhm_arr': fwhm_arr,
                                      'amp_arr': amplitude_arr,
                                      'lock': lock,
                                      'out_dict': cache,
                                      'tag': tag})
            p.start()
            process_list.append(p)
            if len(process_list) >= args.n_proc:
                for p in process_list:
                    p.join()
                print('did batch in %e' % ((time.time()-t_start)/60.0))
                process_list = []

for p in process_list:
    p.join()

with open(args.out_file, 'wb') as file_handle:
    np.savez(file_handle, **cache)

print('that took %e' % (time.time()-t_start))
