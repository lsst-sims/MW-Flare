from __future__ import with_statement

import numpy as np
import os

from mdwarf_utils import light_curve_from_params, generate_light_curve_params
from mdwarf_utils import light_curve_from_class
import time

import multiprocessing as mproc

import argparse

def get_clean_dexes(t_arr, f_arr, t_peak_arr, fwhm_min_arr):
    # do something more rigorous with np.interp
    fixed_dexes = []
    fwhm_arr = fwhm_min_arr/(24.0*60.0)
    fixed_dexes.append(0)
    fixed_dexes.append(len(t_arr)-1)
    for i_peak in range(len(t_peak_arr)):
        peak_dex = np.argmin(np.abs(t_arr-t_peak_arr[i_peak]))
        beginning_dex = np.argmin(np.abs(t_arr-t_peak_arr[i_peak]+fwhm_arr[i_peak]))
        ending_dex = np.argmin(np.abs(t_arr-t_peak_arr[i_peak]+9.0*fwhm_arr[i_peak]))
        fixed_dexes.append(peak_dex)
        fixed_dexes.append(beginning_dex)
        fixed_dexes.append(ending_dex)

    median_flux = np.median(f_arr)
    print('median %e' % median_flux)

    rtol = 0.01
    mtol = 0.001

    keep_going = True
    out_dexes = np.array(range(0,len(t_arr)))
    while keep_going:
        keep_going = False
        dex_dexes = range(0,len(out_dexes),2)
        omitted_dexes = range(1,len(out_dexes),2)
        omitted_dexes = out_dexes[omitted_dexes]
        trial_dexes = list(out_dexes[dex_dexes]) + fixed_dexes
        trial_dexes_arr = np.sort(np.unique(np.array(trial_dexes)))
        f_interped = np.interp(t_arr, t_arr[trial_dexes_arr], f_arr[trial_dexes_arr])
        d_flux = np.abs(f_arr-f_interped)
        bad_dexes = np.where(d_flux>rtol*f_arr+mtol*median_flux)
        for dex in bad_dexes[0]:
            if dex in omitted_dexes:
                trial_dexes.append(dex)

        trial_dexes = np.sort(np.unique(np.array(trial_dexes)))
        if len(trial_dexes) < len(out_dexes):
            keep_going = True
        out_dexes = trial_dexes
        print('keeping %d of %d -- %d' % (len(out_dexes),len(t_arr),len(t_peak_arr)))

    final_pass = True
    while final_pass:
        f_interped = np.interp(t_arr,t_arr[out_dexes],f_arr[out_dexes])
        d_flux = np.abs(f_arr-f_interped)
        bad_dexes = np.where(d_flux>rtol*f_arr+mtol*median_flux)
        print('bad dexes %d' % len(bad_dexes[0]))
        out_dexes = list(out_dexes)
        if len(bad_dexes[0])==0:
            final_pass = False
        for dex in bad_dexes[0]:
            out_dexes.append(dex)
        out_dexes = np.sort(np.unique(np.array(out_dexes)))

    print('keeping %d of %d -- %d' % (len(out_dexes),len(t_arr),len(t_peak_arr)))

    return out_dexes


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
