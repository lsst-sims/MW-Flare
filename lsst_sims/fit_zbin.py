from __future__ import with_statement
import os
import numpy as np
import copy
import argparse

from gaussian_process import ExpSquaredKernel, Covariogram, GaussianProcess

class AdaptiveExpSquaredKernel(ExpSquaredKernel):

    def __init__(self, training_pts=None):
        """
        dim=1 means there is one length scale hyper parameter
        dim=N means that each dimension has its own length scale
        """
        self._training_pts = training_pts
        super(AdaptiveExpSquaredKernel, self).__init__(dim=1)

    def _eval(self, pt1, pt2_list):
        if isinstance(pt1, float):
            ddsq = np.power((pt1-pt2_list), 2)
            med = np.median(np.power(pt1-self._training_pts,2))
        else:
            ddsq = np.power((pt1-pt2_list),2).sum(axis=1)
            med = np.median(np.power(pt1-self._training_pts,2).sum(axis=1), axis=1)

        if med<1.0e-20:
            print 'med ',med
            print 'pt ',pt1

            print 'ddsq ',ddsq
            exit(1)

        return np.exp(-0.5*ddsq/(med*np.power(self._hyper_params,2)))


class UnityGP(GaussianProcess):
    def mean_fn(self, pt):
        if not hasattr(self, '_fn_med'):
            self._fn_med = np.median(self.training_fn)
        return self._fn_med

class LastThreeMeanGP(GaussianProcess):

    def _calc_mean(self, pt):
        if not hasattr(self, '_sqrt_nugget'):
            if len(self.covariogram.nugget)!=len(self.training_pts):
                raise RuntimeError("have not assigned nugget yet %d %d" %
                                   (len(self.training_pts), len(self.nugget)))
            self._sqrt_nugget = 1.0/np.sqrt(self.covariogram.nugget)

        dd = np.abs(pt-self.training_pts)
        sorted_dex = np.argsort(dd)
        if dd[sorted_dex[0]]>0.0001:
            wanted_dex = sorted_dex[:3]
        else:
            wanted_dex = sorted_dex[1:4]

        mean_num = (self.training_fn[wanted_dex]*self._sqrt_nugget[wanted_dex]).sum()
        mean_denom = self._sqrt_nugget[wanted_dex].sum()
        return mean_num/mean_denom

    def mean_fn(self, pt_list):

        if isinstance(pt_list, float):
            return self._calc_mean(pt_list)
        else:
            output = []
            for pp in pt_list:
                output.append(self._calc_mean(pp))
            return np.array(output)


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
        gp = UnityGP(covar)
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
