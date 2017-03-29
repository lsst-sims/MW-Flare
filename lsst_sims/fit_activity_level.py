from __future__ import with_statement

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse
import numpy as np
import os

__all__ = ["find_fraction_spec_active"]

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


def _build_activity_model():
    """
    Read in data taken from West et al. 2008 (AJ 135, 785).
    For each spectral type, return the parametrs A, tau, B
    needed to model spectral activity as a function of
    distance from the Galactic Plane as

    fraction active = A*exp(-z/tau) + B
    """
    data_dir = 'data/activity_by_type'

    dtype = np.dtype([('z', float), ('frac', float),
                      ('min', float), ('max', float)])

    model_aa = []
    model_bb = []
    model_tau = []
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
        model_aa.append(aa)
        model_bb.append(bb)
        model_tau.append(tau)

    return np.array(model_aa), np.array(model_bb), np.array(model_tau)

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
    that are active
    """

    if not hasattr(find_fraction_spec_active, '_model_aa'):
        (find_fraction_spec_active._model_aa,
         find_fraction_spec_active._model_bb,
         find_fraction_spec_active._model_tau) = _build_activity_model()

    i_star_type = int(star_type[1])

    aa = find_fraction_spec_active._model_aa[i_star_type]
    tau = find_fraction_spec_active._model_tau[i_star_type]
    bb = find_fraction_spec_active._model_bb[i_star_type]
    return aa*np.exp(-1.0*z/tau) + bb


def _find_fudge_factor():
    """
    Find the 'fudge factor' relating the scale height of spectral
    activity to the scale height of flaring activity from Figure 12
    of Hilton et al 2010 (AJ 140, 1402)
    """
    data_dir = "data"
    dtype = np.dtype([('z', float), ('frac', float)])

    active_data = np.genfromtxt(os.path.join(data_dir,
                                             'activity_rate_Hilton_et_al_2010.txt'),
                                dtype=dtype)

    flare_data = np.genfromtxt(os.path.join(data_dir,
                                            'flare_rate_Hilton_et_al_2010.txt'),
                               dtype=dtype)


    tau_grid = np.arange(0.1, 200.0, 0.1)
    offset_grid = np.arange(1.0, 200.0)

    error_best_flare = None
    error_best_active = None
    for tau in tau_grid:
        for offset in offset_grid:
            active_model = 1.0-np.exp(-1.0*(active_data['z']-offset)/tau)
            active_error = np.power(active_model-active_data['frac'],2).sum()
            if error_best_active is None or active_error<error_best_active:
                error_best_active = active_error
                tau_active = tau
                offset_active = offset

            flare_model = 1.0-np.exp(-1.0*(flare_data['z']-offset)/tau)
            flare_error = np.power(flare_model-flare_data['frac'],2).sum()
            if error_best_flare is None or flare_error<error_best_flare:
                error_best_flare = flare_error
                tau_flare = tau
                offset_flare = offset


    print 'tau_active: %.9g; %.2e' % (tau_active, error_best_active)
    print 'tau_flare: %.9g; %.2e' % (tau_flare, error_best_flare)
    print 'tau_flare/tau_active: %.9g' % (tau_flare/tau_active)
    print 'offset_active: %.9g' % offset_active
    print 'offset_flare: %.9g' % offset_flare
    return tau_flare/tau_active, tau_flare, offset_flare, tau_active, offset_active

def find_fraction_flare_active(star_type, z):
    """
    Find the fraction of a spectral type that is flaring active
    by finding a model for spectral activity based on the data
    in West et al. 2008 (AJ 135, 785) and scaling the scale
    height of activity by a 'fudge factor' determined by comparing
    the cumulative distributions of spectrally and flaring active
    stars in Figure 12 of Hilton et al 2010 (AJ 140, 1402)

    Parameters
    ----------
    star_type is a string indicating the star's spectral type (M0-M8)

    z is the star's distance from the Galactic Plane in parsecs

    Returns
    -------
    The fraction of that spectral type at that Galactic Plane distance
    that are active
    """
    if not hasattr(find_fraction_flare_active, '_spec_model_aa'):
        (find_fraction_flare_active._spec_model_aa,
         find_fraction_flare_active._spec_model_bb,
         find_fraction_flare_active._spec_model_tau) = _build_activity_model()
        params = _find_fudge_factor()
        find_fraction_flare_active._fudge_factor = params[0]

    i_star_type = int(star_type[1])

    aa = find_fraction_flare_active._spec_model_aa[i_star_type]
    bb = find_fraction_flare_active._spec_model_bb[i_star_type]
    tau = find_fraction_flare_active._spec_model_tau[i_star_type]
    return aa*np.exp(-1.0*z*tau*find_fraction_flare_active._fudge_factor) + bb



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
