"""
This script will read in the 'fraction of each spectral type that is variable
in each bin in distance from the Galactic plane' data stored in z_bins/ and
try to fit that data with a smooth curve using a Gaussian Process.
"""

import numpy as np
import copy


__all__ = ["KernelBase", "ExpSquaredKernel",
           "Covariogram", "GaussianProcess"]


class KernelBase(object):

    @property
    def n_hyper_params(self):
        return self._n_hyper_params

    @property
    def hyper_params(self):
        return self._hyper_params

    @hyper_params.setter
    def hyper_params(self, val):
        if hasattr(val, '__len__'):
            n_val = len(val)
        else:
            n_val = 1

        if n_val != self._n_hyper_params:
            raise RuntimeError("Trying to give "
                               + "Kernel %d hyper params; " % n_val
                               + "it expects %d" % self._n_hyper_params)

        self._hyper_params = copy.deepcopy(val)

    def __call__(self, pt1, pt2_list):
        if isinstance(pt1, float) and len(pt2_list.shape) != 1:
            raise RuntimeError("passed kernel a float "
                               + "and a %d numpy array" % len(pt2_list.shape))

        if isinstance(pt1, np.ndarray) and len(pt1.shape) != 1:
            raise RuntimeError("passed kernel a %s-D pt" % pt1.shape)

        if isinstance(pt1, np.ndarray) and pt2_list.shape[1] != len(pt1):
            raise RuntimeError("passed kernel a %s-D point " % pt1.shape
                               + " and a %s-D list of points" % pt2_list.shape)

        return self._eval(pt1, pt2_list)


class ExpSquaredKernel(KernelBase):

    def __init__(self, dim=1):
        """
        dim=1 means there is one length scale hyper parameter
        dim=N means that each dimension has its own length scale
        """
        self._n_hyper_params = dim
        self._hyper_params = np.ones(dim, dtype=np.float)

    def _eval(self, pt1, pt2_list):
        """
        pt1 is a single point

        pt2 is a list of points
        """
        if isinstance(pt1, float):
            ddsq = np.power((pt1-pt2_list)/self._hyper_params,2)
        else:
            ddsq = np.power((pt1-pt2_list)/self._hyper_params,2).sum(axis=1)

        return np.exp(-0.5*ddsq)


class CovariogramBase(object):

    @property
    def kriging_param(self):
        return self._kriging_param

    @property
    def kernel(self):
        return self._kernel

    @property
    def covar(self):
        return self._covar

    @property
    def det_covar(self):
        return self._det_covar

    @property
    def covar_inv(self):
        return self._covar_inv

    @property
    def nugget(self):
        return self._nugget

    @nugget.setter
    def nugget(self, val):
        self._nugget = copy.deepcopy(val)

    @property
    def hyper_params(self):
        return np.append(self.kernel.hyper_params, self._kriging_param)

    @hyper_params.setter
    def hyper_params(self, val):
        if len(val) != self._n_hyper_params:
            raise RuntimeError("Trying to give "
                               + "covariogram %d hyper params; " % len(val)
                               + "it expects %d" % self._n_hyper_params)

        self._kernel.hyper_params = copy.deepcopy(val[:-1])
        self._kriging_param = copy.deepcopy(val[-1])

    @property
    def n_hyper_params(self):
        return self._n_hyper_params


class Covariogram(CovariogramBase):

    def __init__(self, kernel):
        self._n_hyper_params = kernel.n_hyper_params + 1
        self._kriging_param = 1.0
        self._kernel = kernel
        self._nugget = 1.0e-5
        self._covar = None
        self._det_covar = None
        self._covar_inv = None
        self._diag = None

    def __call__(self, pt1, pt2_list):
        return self.kernel(pt1, pt2_list)/self.kriging_param

    def build_covar(self, pts_in):
        """
        pts_in should be a 2-D numpy array.  Each row is one of the points
        in parameter space on which we are building the Gaussian Process.

        If we are just constructing a 1-D Gaussian Process, this can be
        a 1-D numpy array in which each point is the abscissa of our
        data points.
        """
        if not isinstance(pts_in, np.ndarray):
            raise RuntimeError("Need to pass a numpy array to Covariogram()")

        self._covar = np.zeros((len(pts_in), len(pts_in)))
        for ix, pt in enumerate(pts_in):
            other_dexes = range(ix, len(pts_in))
            other_pts = pts_in[other_dexes]
            kernel_vals = self(pt, other_pts)
            if self._diag is None:
                if isinstance(self.nugget, float):
                    kernel_vals[0] += self.nugget
                else:
                    kernel_vals[0] += self.nugget[ix]
            else:
                if isinstance(self._diag, float):
                    kernel_vals[0] = self._diag
                else:
                    kernel_vals[0] = self._diag[ix]

            for ii, iy in enumerate(other_dexes):
                self._covar[ix][iy] = kernel_vals[ii]
                if ix != iy:
                    self._covar[iy][ix] = kernel_vals[ii]

        self._det_covar = np.linalg.det(self._covar)

    def assign_diagonal_covar(self, val):
        self._diag = val

    def build_covar_inv(self):
        if self._covar is None:
            raise RuntimeError("Cannot build covar_inv; covar is None")

        self._covar_inv = np.linalg.inv(self._covar)


class GaussianProcess(object):

    def __init__(self, covariogram):
        self._covariogram = covariogram
        self._mean_fn = None
        self._inv_dot_fn = None
        self._training_pts = None
        self._training_fn = None
        self._ln_det = None

    @property
    def covariogram(self):
        return self._covariogram

    @property
    def training_pts(self):
        return self._training_pts

    @property
    def training_fn(self):
        return self._training_fn

    @property
    def mean_fn(self):
        return self._mean_fn

    def build(self, training_pts, training_fn):
        self.covariogram.build_covar(training_pts)
        self.covariogram.build_covar_inv()

        self._mean_fn = np.mean(training_fn)

        self._inv_dot_fn = np.dot(self.covariogram.covar_inv, training_fn-self._mean_fn)

        self._training_pts = copy.deepcopy(training_pts)
        self._training_fn = copy.deepcopy(training_fn)
        self._ln_det = np.log(np.abs(self.covariogram.det_covar))

    def ln_likelihood(self):
        if self._mean_fn is None:
            raise RuntimeError("must call GaussianProcess.build() "
                               "before GaussianProcess.likelihood")

        arg = np.dot(self.training_fn-self._mean_fn, self._inv_dot_fn)
        if arg<0.0:
            return -1.0e20
        return -0.5*arg - 0.5*self._ln_det

    def regress(self, test_pts):
        if self._mean_fn is None:
            raise RuntimeError("must call GaussianProcess.build() "
                               "before GaussianProcess.project()")

        output = []
        for pt in test_pts:
            covar_vals = self.covariogram(pt, self.training_pts)
            ans = np.dot(covar_vals, self._inv_dot_fn)
            output.append(self.mean_fn + ans)

        return np.array(output)
