import unittest
import numpy as np

from gaussian_process import ExpSquaredKernel, Covariogram

class ExpSquaredKernelTest(unittest.TestCase):
    """
    Test that the class implementing the squared exponential
    kernel works properly
    """

    def test_on_scalars(self):
        """
        Test on data that is just scalars
        """
        kernel = ExpSquaredKernel()
        self.assertTrue(kernel.n_hyper_params, 1)
        self.assertEqual(kernel.hyper_params[0], 1.0)
        ans = kernel(8.2, np.array([9.0]))
        control = np.exp(-0.5*np.array([0.8*0.8]))
        np.testing.assert_array_almost_equal(ans/control, np.ones(1), decimal=10)

        pt2_list = np.array([1.3, 1.2, 1.4])
        ans = kernel(1.5, pt2_list)
        control = np.exp(-0.5*np.array([0.04, 0.09, 0.01]))
        np.testing.assert_array_almost_equal(ans/control, np.ones(3), decimal=10)

        kernel.hyper_params = 0.53
        ans = kernel(8.3, np.array([5.4]))
        control = np.exp(-0.5*np.array([np.power(2.9/0.53,2)]))
        np.testing.assert_array_almost_equal(ans/control, np.ones(1), decimal=10)

        ans = kernel(2.3, np.array([1.1, 2.5, 3.7]))
        control = np.exp(-0.5*np.power(np.array([1.2, 0.2, 1.4])/0.53,2))
        np.testing.assert_array_almost_equal(ans/control, np.ones(3), decimal=10)


    def test_on_vectors(self):
        """
        Test on mutli-dimensional points (but only one length scale)
        """
        kernel = ExpSquaredKernel()
        self.assertTrue(kernel.n_hyper_params, 1)
        self.assertEqual(kernel.hyper_params[0], 1.0)

        pt1 = np.array([1.1, 2.4, 5.6])

        pt2_list = np.array([[2.3, 3.1, 5.5]])
        ans = kernel(pt1, pt2_list)
        control = np.array([np.exp(-0.5*(1.2*1.2+0.7*0.7+0.01))])
        np.testing.assert_array_almost_equal(ans/control, np.ones(1), decimal=10)

        pt2_list = np.array([[1.1, 1.1, 1.1], [2.4, 2.4, 2.4]])
        ans = kernel(pt1, pt2_list)
        control = np.exp(-0.5*np.array([1.3*1.3+4.5*4.5, 1.3*1.3+3.2*3.2]))
        np.testing.assert_array_almost_equal(ans/control, np.ones(2), decimal=10)

        kernel.hyper_params = 2.0
        self.assertEqual(kernel.hyper_params, 2.0)

        pt2_list = np.array([[2.3, 3.1, 5.5]])
        ans = kernel(pt1, pt2_list)
        control = np.array([np.exp(-0.5*(1.2*1.2+0.7*0.7+0.01)*0.25)])
        np.testing.assert_array_almost_equal(ans/control, np.ones(1), decimal=10)

        pt2_list = np.array([[1.1, 1.1, 1.1], [2.4, 2.4, 2.4]])
        ans = kernel(pt1, pt2_list)
        control = np.exp(-0.5*np.array([1.3*1.3+4.5*4.5, 1.3*1.3+3.2*3.2])*0.25)
        np.testing.assert_array_almost_equal(ans/control, np.ones(2), decimal=10)

    def test_multi_length_scales(self):
        kernel = ExpSquaredKernel(dim=3)
        kernel.hyper_params = np.array([1.0, 2.0, 3.0])
        np.testing.assert_array_equal(kernel.hyper_params, np.array([1.0, 2.0, 3.0]))

        pt1 = np.array([2.3, 5.6, 8.1])
        pt2_list = np.array([[1.1, 3.4, 7.2]])
        ans = kernel(pt1, pt2_list)
        control = np.exp(-0.5*(1.2*1.2+0.25*2.2*2.2+0.9*0.9/9.0))
        np.testing.assert_array_almost_equal(ans/control, np.ones(1), decimal=10)

        pt2_list = np.array([[1.1, 1.1, 1.1], [2.2, 2.2, 2.2]])
        ans = kernel(pt1, pt2_list)
        control = np.exp(-0.5*np.array([1.2*1.2+0.25*4.5*4.5+7.0*7.0/9.0,
                                        0.1*0.1+3.4*3.4*0.25+5.9*5.9/9.0]))
        np.testing.assert_array_almost_equal(ans/control, np.ones(2), decimal=10)


class CovariogramTestCase(unittest.TestCase):

    def test_nugget_assignment(self):
        kernel = ExpSquaredKernel()
        covariogram = Covariogram(kernel)
        self.assertEqual(covariogram.nugget, 1.0e-5)
        covariogram.nugget=1.0e-4
        self.assertEqual(covariogram.nugget, 1.0e-4)


if __name__ == "__main__":
    unittest.main()
