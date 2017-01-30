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
