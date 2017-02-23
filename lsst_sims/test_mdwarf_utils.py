import unittest
import numpy as np

from mdwarf_utils import prob_of_type
from mdwarf_utils import xyz_from_lon_lat_px

class ProbTestCase(unittest.TestCase):

    longMessage = True

    def test_vectorized(self):
        """
        Test that prob_of_type returns the same values when run on arrays
        of stars as when run on single stars
        """
        rng = np.random.RandomState(813)
        n_stars = 10
        ri = rng.random_sample(n_stars)
        iz = rng.random_sample(n_stars)

        vector_pdf = prob_of_type(ri, iz)
        for ix in range(n_stars):
            single_pdf = prob_of_type(ri[ix], iz[ix])
            for mm in range(8):
                msg = 'failed on star %d; type %d' % (ix, mm)
                self.assertAlmostEqual(vector_pdf[mm][ix]/single_pdf[mm],
                                       1.0, 10, msg=msg)


class XyzTestCase(unittest.TestCase):

    longMessage = True

    def test_vectorize(self):
        """
        Test that vectorized xyz_from_lon_lat_px gives the same results
        as xyz_from_lon_lat_px run one-at-a-time
        """
        rng = np.random.RandomState(813412)
        n_stars = 20
        gal_l = rng.random_sample(n_stars)*360.0
        gal_b = rng.random_sample(n_stars)*180.0 - 90.0
        px = rng.random_sample(n_stars)*1.0 + 0.0001

        vector_xyz = xyz_from_lon_lat_px(gal_l, gal_b, px)

        for ix in range(n_stars):
            xyz = xyz_from_lon_lat_px(gal_l[ix], gal_b[ix], px[ix])
            for ii in range(3):
                msg = 'failed on star %d; dim %d' % (ix, ii)
                self.assertEqual(xyz[ii], vector_xyz[ii][ix], msg=msg)

if __name__ == "__main__":
    unittest.main()
