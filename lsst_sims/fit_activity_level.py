from __future__ import with_statement

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse
import numpy as np
import os

def fit_to_exp_decay(xx_data, yy_data, sigma_data):

    tau_grid = np.arange(xx_data.min(), xx_data.max(), 0.1)
    tau_best = None
    aa_best = None
    bb_best = None
    error_best = None
    sig_term = 1.0/np.power(sigma_data, 2)
    gamma = 1.0/sig_term.sum()
    for tau in tau_grid:
        exp_term = np.exp(-1.0*xx_data/tau)
        theta = (exp_term*sig_term).sum()

        aa_num = (exp_term*yy_data*sig_term).sum() - theta*gamma*(yy_data*sig_term).sum()
        aa_denom = (np.exp(-2.0*xx_data/tau)*sig_term).sum() - theta*theta*gamma
        aa = aa_num/aa_denom

        bb = gamma*((yy_data - aa*exp_term)*sig_term).sum()

        err = np.power((yy_data - aa*np.exp(-1.0*xx_data/tau) - bb)/sigma_data, 2).sum()
        if error_best is None or err<error_best:
            error_best = err
            aa_best = aa
            bb_best = bb
            tau_best = tau

    return aa_best, tau_best, bb_best


def find_fraction_spec_active(star_type, z):
    """
    Find the fraction of a spectral type that is active (in the spectroscopic
    sense of  as a function of West et al. 2008 (AJ 135, 785), at a given
    distance from the Galactic Plane

    Parameters
    ----------
    star_type is a string indicating the star's spectral type (M0-M8)

    z is the star's distance from the Galactic Plane in parsecs

    Returns
    -------
    The fraction of that spectral type at that Galactic Plane distance
    that are activve
    """

    if not hasattr(find_fraction_spec_active, '_model_dict'):
        model_dict = {}
        data_dir = 'data/activity_by_type'

        dtype = np.dtype([('z', float), ('frac', float),
                          ('min', float), ('max', float)])

        for i_type in range(9):
            model_type = 'M%d' % i_type
            data_name = os.path.join(data_dir, '%s.txt' % model_type)
            data = np.genfromtxt(data_name, dtype=dtype)

            sigma_arr = []
            for nn, xx in zip(data['min'], data['max']):
                if nn>1.0e-20 and xx<0.999:
                    sigma = 0.5*(xx-nn)
                else:
                    sigma = xx-nn
                sigma_arr.append(sigma)
            aa, tau, bb = fit_to_exp_decay(data['z'], data['frac'], np.array(sigma_arr))
            model_dict[model_type] = {}
            model_dict[model_type]['a'] = aa
            model_dict[model_type]['tau'] = tau
            model_dict[model_type]['b'] = bb

            find_fraction_spec_active._model_dict = model_dict

    aa = find_fraction_spec_active._model_dict[star_type]['a']
    tau = find_fraction_spec_active._model_dict[star_type]['tau']
    bb = find_fraction_spec_active._model_dict[star_type]['b']
    return aa*np.exp(-1.0*z/tau) + bb


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, default='type_fits')

    args = parser.parse_args()
    if args.outdir is None:
        raise RuntimeError("need to specify an output directory")

    type_list = ['M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8']
    data_dir = 'data/activity_by_type'

    dtype = np.dtype([('z', float), ('frac', float),
                      ('min', float), ('max', float)])

    plot_dir = 'plots'
    plt.figsize = (30, 30)

    xx_test = np.arange(0.0, 1000.0, 1.0)

    for i_fig, spec_type in enumerate(type_list):

        data_name = os.path.join(data_dir, '%s.txt' % spec_type)
        data = np.genfromtxt(data_name, dtype=dtype)

        yy_test = find_fraction_spec_active(spec_type, xx_test)

        plt.subplot(3,3,i_fig+1)
        plt.errorbar(data['z'], data['frac'],
                     yerr = np.array([data['frac']-data['min'],
                                      data['max']-data['frac']]),
                     marker='o', linestyle='')
        plt.plot(xx_test, yy_test, color='r')
        plt.title(spec_type, fontsize=10)
        plt.ylim(0.0, 1.1)
        yticks = np.arange(0.0, 1.15, 0.1)
        ylabels = ['%.1f' % xx if ii%4==0 else ''
                   for ii, xx in enumerate(yticks)]
        plt.yticks(yticks,ylabels)

    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'exp_decay_by_type.png'))
