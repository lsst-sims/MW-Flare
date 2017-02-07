from __future__ import with_statement
import os
import numpy as np
import copy
import argparse

from gaussian_process import ExpSquaredKernel, Covariogram, GaussianProcess

class LinearMeanGP(GaussianProcess):

    def __init__(self,covariogram,pt1,pt2,asymptote):

        self._slope = (pt1[1]-pt2[1])/(pt1[0]-pt2[0])
        self._intercept = pt1[1] - self._slope*pt1[0]
        self._max = pt2[0]
        self._asymptote = asymptote

        super(LinearMeanGP, self).__init__(covariogram)


    def mean_fn(self, pt_list):
        if isinstance(pt_list, float):
            pts = np.array([pt_list])
        else:
            pts = pt_list

        return self._intercept + self._slope*pt_list


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', type=str, default=None)
    parser.add_argument('--outfile', type=str, default=None)
    parser.add_argument('--anchor', type=float, nargs='+', default=None)

    args = parser.parse_args()
    if args.infile is None or args.outfile is None:
        raise RuntimeError("input %s output %s" % (args.infile, args.outfile))

    dtype = np.dtype([('class', float), ('frac', float),
                      ('minus', float), ('plus', float)])

    data = np.genfromtxt(args.infile, dtype=dtype)

    nugget = np.zeros(len(data))
    for ix in range(len(data)):
        if data['minus'][ix]>1.0e-6 and data['plus'][ix]<0.99999:
            nugget[ix] = 0.5*(data['plus'][ix]-data['minus'][ix])
        elif data['minus'][ix]<1.0e-6:
            nugget[ix] = data['plus'][ix]-data['frac'][ix]
        else:
            nugget[ix] = data['frac'][ix]-data['minus'][ix]

    nugget = np.power(nugget, 2)

    kernel = ExpSquaredKernel(dim=1)
    covar = Covariogram(kernel)
    covar.nugget = nugget
    if args.anchor is None:
        gp = GaussianProcess(covar)
    else:
        gp = LinearMeanGP(covar,
                          (args.anchor[0], args.anchor[1]),
                          (args.anchor[2], args.anchor[3]),
                          args.anchor[4])

    min_hyper_params = np.array([0.5, 1.0])
    gp.covariogram.hyper_params = min_hyper_params


    gp.build(data['class'], data['frac'])
    max_like = None
    rng = np.random.RandomState(44)
    for ell in np.arange(0.05, 3.0, 0.05):
        for kk in  np.arange(0.1, 100.0, 0.1):
            test_params = np.array([ell, kk])
            gp.covariogram.hyper_params = test_params
            gp.build(data['class'], data['frac'])
            ln_like = gp.ln_likelihood()
            if max_like is None or ln_like > max_like:
                print ln_like, test_params, -2.0*(ln_like+0.5*gp._ln_det)
                max_like = ln_like
                min_hyper_params = copy.deepcopy(test_params)
            else:
                gp.covariogram.hyper_params = min_hyper_params


    gp.covariogram.hyper_params = min_hyper_params
    gp.build(data['class'], data['frac'])
    print gp.ln_likelihood(), gp._ln_det

    x_test = np.arange(0.0,16.0, 0.1)
    y_test = gp.regress(x_test)

    with open(args.outfile, "w") as output_file:
        output_file.write('# data %s\n' % args.infile)
        output_file.write('# hyper params: %.3f %.3f\n' %
                          (gp.covariogram.hyper_params[0],
                           gp.covariogram.hyper_params[1]))
        for xx, yy in zip(x_test, y_test):
            output_file.write("%e %e\n" % (xx, yy))
