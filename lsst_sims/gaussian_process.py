"""
This script will read in the 'fraction of each spectral type that is variable
in each bin in distance from the Galactic plane' data stored in z_bins/ and
try to fit that data with a smooth curve using a Gaussian Process.
"""

import numpy as np
import copy


__all__ = ["KernelBase", "ExpSquaredKernel", "Covariogram"]


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

        self._hyper_params = val

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
    def covar_inv(self):
        return self._covar_inv

    @property
    def nugget(self):
        return self._nugget

    @nugget.setter
    def nugget(self, val):
        if not isinstance(val, float):
            raise RuntimeError("Nugget should be a float; "
                               "passing %s" % type(val))
        self._nugget = val

    @property
    def hyper_params(self):
        return np.append(self.kernel.hyper_params, self._kriging_param)

    @hyper_params.setter
    def hyper_params(self, val):
        if len(val) != self._n_hyper_params:
            raise RuntimeError("Trying to give "
                               + "covariogram %d hyper params; " % len(val)
                               + "it expects %d" % self._n_hyper_params)

        self._kernel.hyper_params = val[:-1]
        self._kriging_param = val[-1]

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
        self._covar_inv = None

    def __call__(self, pts_in):
        """
        pts_in should be a 2-D numpy array.  Each row is one of the points
        in parameter space on which we are building the Gaussian Process.

        If we are just constructing a 1-D Gaussian Process, this can be
        a 1-D numpy array in which each point is the abscissa of our
        data points.
        """
        if not isinstance(pts_in, np.ndarray):
            raise RuntimeError("Need to pass a numpy array to Covariogram()")

        self._covar = []
        for ix, pt in enumerate(pts_in):
            kernel_vals = self.kernel(pt, pts_in)
            kernel_vals[ix] += self.nugget
            self._covar.append(kernel_vals)

        self._covar = np.array(self._covar)/self.kriging_param
        self._covar_inv = np.linalg.inv(self._covar)
