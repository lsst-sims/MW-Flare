from __future__ import with_statement
import os
import numpy as np
import copy

import argparse


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

    mean_val = np.mean(data['frac'])

    sigma = np.zeros(len(data))
    for ix in range(len(data)):
        if data['minus'][ix]>1.0e-6 and data['plus'][ix]<0.99999:
            sigma[ix] = 0.5*(data['plus'][ix]-data['minus'][ix])
        elif data['minus'][ix]<1.0e-6:
            sigma[ix] = data['plus'][ix]-data['frac'][ix]
            #data['frac'][ix] = 0.5*(data['frac'][ix]+data['plus'][ix])
        else:
            sigma[ix] = data['frac'][ix]-data['minus'][ix]
            #data['frac'][ix] = 0.5*(data['frac'][ix]+data['minus'][ix])

    like_best = None
    for ell in np.arange(0.05, 5.05, 0.05):
        covar = np.zeros((len(data), len(data)))
        for ix in range(len(data)):
            for iy in range(len(data)):
                if ix==iy:
                    covar[ix][iy] = sigma[ix]*sigma[ix]
                else:
                    dd = np.abs(data['class'][ix]-data['class'][iy])
                    covar[ix][iy] = sigma[ix]*sigma[iy]*np.exp(-0.5*np.power(dd/ell,2))

        covar_inv = np.linalg.inv(covar)
        covar_det = np.linalg.det(covar)
        chisq = np.dot(data['frac']-mean_val,
                       np.dot(covar_inv, data['frac']-mean_val))
        ln_like = -0.5*chisq-0.5*np.log(np.abs(covar_det))

        if like_best is None or ln_like>like_best:
            print 'ell %e ln_like %e' % (ell, ln_like)
            like_best = ln_like
            ell_best = ell
            covar_best = covar
            covar_inv_best = covar_inv

    covar = covar_best
    covar_inv = covar_inv_best
    ell = ell_best

    x_test = np.arange(0.0, 16.0, 0.1)
    with open(args.outfile, 'w') as output_file:
        for xx in x_test:
            dd = np.abs(xx-data['class'])
            covar_star = sigma*sigma*np.exp(-0.5*np.power(dd/ell,2))
            yy = mean_val + np.dot(covar_star,
                                   np.dot(covar_inv, data['frac']-mean_val))
            output_file.write('%e %e\n' % (xx,yy))

    print 'mean_val ',mean_val
    chisq = np.dot(data['frac']-mean_val,
                   np.dot(covar_inv, data['frac']-mean_val))
    ln_like = -0.5*chisq-0.5*np.log(np.abs(covar_det))
