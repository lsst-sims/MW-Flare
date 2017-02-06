from __future__ import with_statement
import os
import numpy as np
import copy
import argparse

from gaussian_process import ExpSquaredKernel, Covariogram, GaussianProcess


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', type=str, default=None)
    parser.add_argument('--outfile', type=str, default=None)

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
    gp = GaussianProcess(covar)

    min_hyper_params = np.array([0.5, 1.0])
    gp.covariogram.hyper_params = min_hyper_params


    gp.build(data['class'], data['frac'])
    max_like = gp.ln_likelihood()
    rng = np.random.RandomState(44)
    for ix in range(50000):
        delta = (rng.random_sample(2)-0.5)
        test_params = min_hyper_params + delta
        gp.covariogram.hyper_params = test_params
        gp.build(data['class'], data['frac'])
        ln_like = gp.ln_likelihood()
        if ln_like > max_like:
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

    print gp.covariogram.covar
    print '\n\n'
    print gp.covariogram(0.0,gp.training_pts)
    print '\n\n'
    print gp.training_pts
    print '\n\n'
    print gp._inv_dot_fn

    with open(args.outfile, "w") as output_file:
        for xx, yy in zip(x_test, y_test):
            output_file.write("%e %e\n" % (xx, yy))
