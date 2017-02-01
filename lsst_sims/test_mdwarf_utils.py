import unittest
import numpy as np

from mdwarf_utils import prob_of_type

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

if __name__ == "__main__":
    unittest.main()
