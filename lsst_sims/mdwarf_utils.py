import numpy as np
from lsst.sims.utils import radiansFromArcsec

__all__ = ["xyz_from_lon_lat_px"]

def xyz_from_lon_lat_px(lon, lat, px):
    """
    Parameters
    ----------
    lon is galactic longitude in degrees

    lat is galactic latitude in degrees

    px in parallax in arcseconds

    Returns
    -------
    A numpy array containing the vector from Sgr A* to the
    input star.  Distances in parsecs.
    """


    _d_center = 8.33e3 # distance from the Sun to Sgr A* from
                       # Gillessen et al 2009 (ApJ 692, 1075) eqn 11

    _lon_center = 359.9442
    _lat_center = -0.0462 # galactic latitude and longitude of Sgr A* from
                          # http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=Sagittarius+A&submit=SIMBAD+search

    _au_to_parsec = 1.0/206265.0

    if not hasattr(xyz_from_lon_lat_px, '_xyz_center'):
        lon_rad = np.radians(_lon_center)
        lat_rad = np.radians(_lat_center)
        xx = np.cos(lon_rad)*np.cos(lat_rad)
        yy = np.sin(lon_rad)*np.cos(lat_rad)
        zz = np.sin(lat_rad)

        xyz_from_lon_lat_px._xyz_center = np.array([_d_center*xx,
                                                    _d_center*yy,
                                                    _d_center*zz])

    lon_rad = np.radians(lon)
    lat_rad = np.radians(lat)
    dd = _au_to_parsec/radiansFromArcsec(px)

    xx = dd*np.cos(lon_rad)*np.cos(lat_rad)
    yy = dd*np.sin(lon_rad)*np.cos(lat_rad)
    zz = dd*np.sin(lat_rad)

    return np.array([xx-xyz_from_lon_lat_px._xyz_center[0],
                     yy-xyz_from_lon_lat_px._xyz_center[1],
                     zz-xyz_from_lon_lat_px._xyz_center[2]])


def prob_of_type(r_i, i_z):
    """
    Find the relative probabilities that a star is of spectral types
    M0-M7 as detailed in the caption for Table 1 of Kowalski et al 2009
    (ApJ 138, 633)

    Parameters
    ----------
    r_i is the star's SDSS r-i color

    i_z is the star's SDSS i-z color

    Returns
    -------
    A numpy array containing the probability density value of the star's color
    in the 2-D Gaussian PDF associated with each spectral subtype (M0-M7).
    If more than one star was passed in, then each row of the numpy array
    will corresponds to a different stellar type, i.e.

    output[1][2]

    will be the probability that star 2 is of type M1
    """

    if not hasattr(prob_of_type, 'r_i'):
        dtype = np.dtype([('type', str, 2),
                          ('r_i', float), ('i_z', float),
                          ('cov_00', float), ('cov_01', float),
                          ('cov_10', float), ('cov_11', float)])

        input_data = np.genfromtxt('color_covar_data.txt', dtype=dtype)
        prob_of_type.r_i = []
        prob_of_type.i_z = []
        prob_of_type.covar_inv = []
        prob_of_type.sqrt_det = []
        for ix, row in enumerate(input_data):
            assert row[0] == 'M%d' % ix
            prob_of_type.r_i.append(row[1])
            prob_of_type.i_z.append(row[2])

            covar = np.array([[row[3], row[4]], [row[5], row[6]]])
            covar_inv = np.linalg.inv(covar)
            prob_of_type.covar_inv.append(covar_inv)
            sqrt_det = np.sqrt(np.linalg.det(covar))
            prob_of_type.sqrt_det.append(sqrt_det)


    output = []
    for ix in range(len(prob_of_type.r_i)):
        x_minus_mu = np.array([r_i-prob_of_type.r_i[ix],
                               i_z-prob_of_type.i_z[ix]]).transpose()

        if len(x_minus_mu.shape)>1:
            arg1 = np.dot(x_minus_mu, prob_of_type.covar_inv[ix])

            # see http://stackoverflow.com/questions/15616742/vectorized-way-of-calculating-row-wise-dot-product-two-matrices-with-scipy
            arg = np.einsum('ij,ij->i',arg1,x_minus_mu)
        else:
            arg = np.dot(x_minus_mu, np.dot(prob_of_type.covar_inv[ix], x_minus_mu))

        exp_term = np.exp(-0.5*arg)
        output.append(exp_term/(2.0*np.pi*prob_of_type.sqrt_det[ix]))

    return np.array(output)
