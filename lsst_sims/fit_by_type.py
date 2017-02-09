from __future__ import with_statement

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse
import numpy as np
import os

from gaussian_process import ExpSquaredKernel, Covariogram, GaussianProcess

class LinearMeanGP(GaussianProcess):

    def mean_fn(self, pt_list):
        if not hasattr(self, '_slope'):
            self._slope = (self.training_fn[-1]-self.training_fn[0])
            self._slope /= (self.training_pts[-1]-self.training_pts[0])

            self._intercept = self.training_fn[0] - self._slope*self.training_pts[0]

        return self._slope*pt_list + self._intercept

class ExpMeanGP(GaussianProcess):

    def _fit_tau(self):
        floor_grid = np.arange(0.01, 0.3, 0.01)
        tau_grid = np.arange(10.0, 10000.0, 10.0)
        norm_grid = np.arange(0.0, 1.0, 0.05)
        chisq_best = None
        for tau in tau_grid:
            for norm in norm_grid:
                for floor in floor_grid:
                    yy = norm*np.exp(-1.0*self.training_pts/tau) + floor
                    chisq = (np.power(yy-self.training_fn,2)/self.covariogram.nugget).sum()
                    if chisq_best is None or chisq<chisq_best:
                        chisq_best = chisq
                        self._tau = tau
                        self._norm = norm
                        self._floor = floor
        print self._norm, self._tau, self._floor


    def mean_fn(self, pt_list):
        if not hasattr(self, '_tau'):
            self._fit_tau()

        return self._norm*np.exp(-1.0*pt_list/self._tau) + self._floor

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, default=None)

    args = parser.parse_args()
    if args.outdir is None:
        raise RuntimeError("need to specify an output directory")

    type_list = ['M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8']
    data_dir = 'data/activity_by_type'

    dtype = np.dtype([('z', float), ('frac', float),
                      ('min', float), ('max', float)])

    plot_dir = 'plots'
    plt.figsize = (30, 30)

    for i_fig, spec_type in enumerate(type_list):
        data_name = os.path.join(data_dir, '%s.txt' % spec_type)
        data = np.genfromtxt(data_name, dtype=dtype)

        nugget = []
        for nn, xx in zip(data['min'], data['max']):
            if nn>1.0e-20 and xx<0.999:
                sigma = 0.5*(xx-nn)
            else:
                sigma = xx-nn
            nugget.append(np.power(sigma,2))

        nugget = np.array(nugget)

        kernel = ExpSquaredKernel(dim=1)
        covariogram = Covariogram(kernel)
        covariogram.nugget = nugget

        gp = ExpMeanGP(covariogram)
        max_like = None
        for ell in np.arange(10.0, 700.0, 10.0):
            for log_krig in np.arange(-3.0, 6.1, 0.1):
                test_params = np.array([ell, np.power(10, log_krig)])
                gp.covariogram.hyper_params = test_params
                gp.build(data['z'], data['frac'])
                like = gp.ln_likelihood()
                if max_like is None or like>max_like:
                    max_like = like
                    best_params = test_params

        print spec_type,best_params,max_like
        gp.covariogram.hyper_params = best_params
        gp.build(data['z'], data['frac'])

        out_name = os.path.join(args.outdir, '%s_gp.txt' % spec_type)
        with open(out_name, 'w') as out_file:
            xx_test = np.arange(0.0, 1000.0, 10.0)
            yy_test = gp.regress(xx_test)
            for xx, yy in zip(xx_test, yy_test):
                out_file.write('%e %e\n' % (xx, yy))

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
    plt.savefig(os.path.join(plot_dir, 'gp_by_type.png'))
