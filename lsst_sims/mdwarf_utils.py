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
    A dict containing the probability density value of the star's color
    in the 2-D Gaussian PDF associated with each spectral subtype (M0-M7).
    """

    if not hasattr(prob_of_type, 'data'):
        dtype = np.dtype([('type', str, 2),
                          ('r_i', float), ('i_z', float),
                          ('cov_00', float), ('cov_01', float),
                          ('cov_10', float), ('cov_11', float)])

        input_data = np.genfromtxt('color_covar_data.txt', dtype=dtype)
        prob_of_type.data = {}
        for row in input_data:
            prob_of_type.data[row[0]] = {}
            prob_of_type.data[row[0]]['r_i'] = row[1]
            prob_of_type.data[row[0]]['i_z'] = row[2]

            covar = np.array([[row[3], row[4]], [row[5], row[6]]])
            covar_inv = np.linalg.inv(covar)
            prob_of_type.data[row[0]]['covar_inv'] = covar_inv
            sqrt_det = np.sqrt(np.linalg.det(covar))
            prob_of_type.data[row[0]]['sqrt_det'] = sqrt_det


    output = {}
    for name in prob_of_type.data:
        x_minus_mu = np.array([r_i-prob_of_type.data[name]['r_i'],
                               i_z-prob_of_type.data[name]['i_z']])

        arg = np.dot(x_minus_mu,
                     np.dot(prob_of_type.data[name]['covar_inv'], x_minus_mu))

        exp_term = np.exp(-0.5*arg)
        output[name] = exp_term/(2.0*np.pi*prob_of_type.data[name]['sqrt_det'])

    return output
