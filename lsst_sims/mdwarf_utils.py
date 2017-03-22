from __future__ import with_statement
import numpy as np
import os
from lsst.sims.utils import radiansFromArcsec
from lsst.sims.photUtils import Bandpass, BandpassDict
from lsst.utils import getPackageDir

__all__ = ["xyz_from_lon_lat_px", "prob_of_type", "draw_energies",
           "duration_from_energy", "fwhm_from_duration",
           "amplitude_from_fwhm_energy", "lsst_flare_fluxes_from_u"]

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
    M0-M9 as detailed in the caption for Table 2 of West et al 2011
    (AJ 141, 97)

    Parameters
    ----------
    r_i is the star's SDSS r-i color

    i_z is the star's SDSS i-z color

    Returns
    -------
    A numpy array containing the probability density value of the star's color
    in the 2-D Gaussian PDF associated with each spectral subtype (M0-M9).
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

        input_data = np.genfromtxt('data/color_covar_data_West_et_al_2011.txt',
                                   dtype=dtype)
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


def draw_energies(stellar_class, duration, rng):
    """
    Parameters
    ----------
    stellar_class is a string string denoting the class of flare star.  Must be
    selected from
        'early_inactive'
        'early_active'
        'mid_inactive'
        'mid_active'
        'late_active'

    duration is a float indicating the length of time in days that you want to
    simulate.

    rng is an instantiation of numpy.random.RandomState that will be used as
    a random number generator

    Returns
    -------
    A numpy array of times at which flares occurred (in days)
    A numpy array of the energies of the flares (in ergs in the Johnson U band)
    """

    duration_hours = 24.0*duration

    if not hasattr(draw_energies, '_params_dict'):
        draw_energies._params_dict = {}

        param_file = os.path.join('lsst_sims', 'data', 'hilton_phd_table_4.3.txt')
        if not os.path.exists(param_file):
            param_file = os.path.join('data', 'hilton_phd_table_4.3.txt')
            if not os.path.exists(param_file):
                raise RuntimeError("Cannot find hilton_phd_table_4.3.txt.\n"
                                   "Make sure you are calling draw_energies from "
                                   "MW-Flare or MW-Flare/lsst_sims")

        with open(param_file, 'r') as input_file:
            input_lines = input_file.readlines()
            for line in input_lines:
                if line[0] != '#':
                    vv = line.split()
                    _class = vv[0]
                    draw_energies._params_dict[_class] = {}
                    draw_energies._params_dict[_class]['alpha'] =  float(vv[1])
                    draw_energies._params_dict[_class]['beta'] = float(vv[2])
                    draw_energies._params_dict[_class]['log(emin)'] = float(vv[3])
                    draw_energies._params_dict[_class]['log(emax)'] = float(vv[4])

    alpha = draw_energies._params_dict[stellar_class]['alpha']
    beta = draw_energies._params_dict[stellar_class]['beta']
    logemin = draw_energies._params_dict[stellar_class]['log(emin)']

    emin = np.power(10.0, logemin)

    # No obvious justification for this value;
    # it just seems reasonable
    e_abs_max = np.power(10.0,34)

    total_per_hour = np.power(10.0, alpha+beta*logemin)
    total_n_flares = int(np.round(total_per_hour*duration_hours))

    uniform_deviate = rng.random_sample(total_n_flares)

    energy_list = emin*np.power(1.0-uniform_deviate, 1.0/beta)

    # set any flares that were randomly assigned energy > e_abs_max
    # to == e_abs_max to avoid randomly getting flares with absurd
    # energies
    energy_list = np.where(energy_list<e_abs_max, energy_list, e_abs_max)

    time_list = rng.random_sample(total_n_flares)*duration_hours

    return time_list/24.0, energy_list

def _set_up_duration_models():
    """
    Compute the mean and standard deviation of log10(duration) in minutes
    in bins of 1 dex in log10(E_kp).

    Data taken from http://github.com/jradavenport/GJ1243-Flares

    Returns
    -------
    numpy array of energies (middle of bins)

    numpy array of mean log10(duration)s

    numpy array of standard deviation of log10(duration)s
    """

    dtype = np.dtype([('start_dex', float), ('stop_dex', float),
                      ('start', float), ('stop', float), ('peak', float),
                      ('rise', float), ('decay', float), ('amp', float),
                      ('e_dur', float), ('e_dur_rise', float),
                      ('e_dur_decay', float), ('flag', float),
                      ('ppl_flare', float), ('ppl_month', float),
                      ('components', float)])

    data = np.genfromtxt('data/gj1243_master_flares.tbl', dtype=dtype)

    duration = (data['stop'] - data['start'])*24.0*60.0

    log_ekp_quiescent = 30.67
    log_ekp = log_ekp_quiescent + np.log10(data['e_dur'])

    e_min_arr = np.array(range(29,33))
    e_max_arr = e_min_arr + 1

    mean_arr = []
    std_dev_arr = []

    for e_min, e_max in zip(e_min_arr, e_max_arr):
        valid = np.where(np.logical_and(log_ekp>=e_min, log_ekp<=e_max))
        bin_data = np.log10(duration[valid])
        mean = np.mean(bin_data)
        std_dev = np.std(bin_data)
        mean_arr.append(mean)
        std_dev_arr.append(std_dev)

    return 0.5*(e_min_arr+e_max_arr), np.array(mean_arr), np.array(std_dev_arr)


def duration_from_energy(energy_u, rng):
    """
    Find duration as a function of energy by drawing from a Gaussian
    whose mean is calculated by linearly fitting the mean duration in
    energy bins from the energy versus duration bin of Figure 10 of
    Hawley et al 2014 (ApJ 797, 121) and whose standard deviation is
    set to 0.16 (also from the data underlying that plot).

    Parameters
    ----------
    energy_u is the energy of the flare in the Johnson U band
    (in ergs)

    rng is an instantiation of numpy.random.RandomState from which
    we will draw our random numbers

    Returns
    -------
    The duration of the flare in minutes, drawn from Figure 10 of
    Hawley et al. 2014 (ApJ 797, 121)
    """

    if not hasattr(duration_from_energy, '_mean_mm'):

        log_ekp, mean, stdev = _set_up_duration_models()

        # transform to Johnson U band
        # See Hawley et al 2014, last
        # paragraph before Section 3
        log_eu = log_ekp + np.log10(0.65)

        duration_from_energy._stdev = 0.16

        n_data = float(len(log_eu))
        mm = (log_eu*(mean-(mean.sum()/n_data))).sum()
        mm = mm/((log_eu*log_eu).sum() - log_eu.sum()*log_eu.sum()/n_data)

        bb = (mean-mm*log_eu).sum()/n_data

        duration_from_energy._mean_mm = mm
        duration_from_energy._mean_bb = bb

    log_e = np.log10(energy_u)

    mean_dur = duration_from_energy._mean_mm*log_e + duration_from_energy._mean_bb
    sigma = duration_from_energy._stdev

    normal_deviate = rng.normal(loc=0.0, scale=1.0, size=len(energy_u))

    duration = normal_deviate*sigma + mean_dur
    return np.power(10.0, duration)


def _f_rise_of_t(t):
    """
    Return the rising flux of a flare as a funciton of time.
    See equation (1) of Davenport et al 2014 (ApJ 797, 122)

    Parameters
    ----------
    t is the time (in units of the fwhm time)

    Returns
    -------
    Flux (in relative units) at those times
    """

    return 1.0 + 1.941*t - 0.175*t*t - 2.246*t*t*t - 1.125*t*t*t*t

def _f_decay_of_t(t):
    """
    Return the decaying flux of a flare as a funciton of time.
    See equation (4) of Davenport et al 2014 (ApJ 797, 122)

    Parameters
    ----------
    t is the time (in units of the fwhm time)

    Returns
    -------
    Flux (in relative units) at those times
    """

    return 0.689*np.exp(-1.6*t) + 0.0303*np.exp(-0.2783*t)

def fwhm_from_duration(dur):
    """
    Parameters
    ----------
    dur is the duration of the flare in minutes.

    Returns
    -------
    The Full Width at Half Maximum time of the flare
    in minutes (this is the independent parameter in
    equations 1 and 4 of Davenport et al 2014
    (ApJ 797, 122).
    """
    return dur/3.0

def amplitude_from_fwhm_energy(t_fwhm, energy_u):
    """
    Based on the analytical flare model of Davenport et al 2014
    (ApJ 797, 122)

    Parameters
    ----------
    t_fwhm is the width in time of the flare when its amplitude is
    one half its maximum (this is the independent parameter in
    equations 1 and 4 of Davenport et al 2014).  In minutes.

    energy is the energy of the flare in ergs in the Johnson U band

    Returns
    -------
    The amplitude (in ergs/s; *not* the relative amplitude) of the flare.
    """

    if not hasattr(amplitude_from_fwhm_energy, '_t_rise'):
        dt = 0.01
        t_rise = np.arange(-1.0, 0.0+0.5*dt, dt)
        t_decay = np.arange(0.0, 7.0, dt)
        f_rise = _f_rise_of_t(t_rise)
        f_decay = _f_decay_of_t(t_decay)
        amplitude_from_fwhm_energy._f_rise = dt*0.5*(f_rise[1:]+f_rise[:-1]).sum()
        amplitude_from_fwhm_energy._f_decay = dt*0.5*(f_decay[1:]+f_decay[:-1]).sum()

    rising_flux = t_fwhm*amplitude_from_fwhm_energy._f_rise
    decaying_flux = t_fwhm*amplitude_from_fwhm_energy._f_decay

    amplitude = energy_u/(60.0*(rising_flux + decaying_flux))

    assert len(amplitude) == len(t_fwhm)
    assert len(amplitude) == len(energy_u)

    return amplitude

def lsst_flare_fluxes_from_u(u_flux):
    """
    Convert from Johnson U band flux to flux in the LSST bands
    by assuming the flare is a 9000K black body (see Section 4
    of Hawley et al 2003, ApJ 597, 535)

    Parameters
    ----------
    flux in Johnson U band (either a float or a numpy array)

    Returns
    -------
    floats/numpy arrays of fluxes in all 6 LSST bands
    """

    if not hasattr(lsst_flare_fluxes_from_u, 'johnson_u_flux'):
        throughputs_dir = getPackageDir('throughputs')
        johnson_dir = os.path.join(throughputs_dir, 'johnson')
        johnson_u_hw = Bandpass()
        johnson_u_hw.readThroughput(os.path.join(johnson_dir, 'johnson_U.dat'))
        atm = Bandpass()
        atm.readThroughput(os.path.join(throughputs_dir,
                                        'baseline', 'atmos_std.dat'))

        wv, sb = johnson_u_hw.multiplyThroughputs(atm.wavelen, atm.sb)
        johnson_u = Bandpass(wavelen=wv, sb=sb)

        boltzmann_k = 1.3807e-16  # erg/K
        planck_h = 6.6261e-27  # erg*s
        _c = 2.9979e10  # cm/s

        hc_over_k = 1.4387e7  # nm*K

        temp = 9000.0  # black body temperature in Kelvin

        bb_wavelen = np.arange(200.0, 1500.0, 0.1)  # in nanometers

        exp_arg = -1.0*hc_over_k/(temp*bb_wavelen)

        log_bb_flambda = -5.0*np.log(bb_wavelen) + exp_arg
        bb_flambda = np.exp(log_bb_flambda)

        # Note: we are ignoring the 2*h*c^2 term, as well as all other
        # normalization terms, because we are only interested in
        # how the fluxes scale relatively between the 7 bands (Johnson U
        # and the 6 LSST bands)

        bb_int_flambda = np.interp(johnson_u.wavelen, bb_wavelen, bb_flambda)

        # calculate the flux (with an arbitrary normalization) of the 9000K
        # black body in the Johnson_U band and the lsst bands

        # integrate the flux of the black body over the Johnson U band
        bb_johnson_u = (0.5*(johnson_u.wavelen[1]-johnson_u.wavelen[0])
                        *(johnson_u.sb[1:]*bb_int_flambda[1:] +
                          johnson_u.sb[:-1]*bb_int_flambda[:-1]).sum())

        lsst_flare_fluxes_from_u.johnson_u_raw_flux = bb_johnson_u

        lsst_flare_fluxes_from_u.lsst_raw_flux_dict = {}
        lsst_bands = BandpassDict.loadTotalBandpassesFromFiles()
        for band_name in ('u', 'g', 'r', 'i', 'z', 'y'):
            bp = lsst_bands[band_name]
            bb_int_flambda = np.interp(bp.wavelen, bb_wavelen, bb_flambda)

            flux = (0.5*(bp.wavelen[1]-bp.wavelen[0])
                    *(bp.sb[1:]*bb_int_flambda[1:] +
                      bp.sb[:-1]*bb_int_flambda[:-1]).sum())

            lsst_flare_fluxes_from_u.lsst_raw_flux_dict[band_name] = flux

    factor = u_flux/lsst_flare_fluxes_from_u.johnson_u_raw_flux
    u_flux = factor*lsst_flare_fluxes_from_u.lsst_raw_flux_dict['u']
    g_flux = factor*lsst_flare_fluxes_from_u.lsst_raw_flux_dict['g']
    r_flux = factor*lsst_flare_fluxes_from_u.lsst_raw_flux_dict['r']
    i_flux = factor*lsst_flare_fluxes_from_u.lsst_raw_flux_dict['i']
    z_flux = factor*lsst_flare_fluxes_from_u.lsst_raw_flux_dict['z']
    y_flux = factor*lsst_flare_fluxes_from_u.lsst_raw_flux_dict['y']

    return (u_flux, g_flux, r_flux, i_flux, z_flux, y_flux)



def light_curve_from_class(stellar_class, years, rng, dt=15.0):
    """
    Simulate a flaring star light curve

    Parameters
    ----------
    stellar_class is either
        early_inactive
        early_active
        mid_inactive
        mid_active
        late_active

    years is the number of years to simulate

    rng is an instantiation of numpy.random.RandomState to be used
    as a random number generator

    dt is the time step to take in seconds

    Returns
    -------
    time is a numpy array containing the times (in days)
    corresponding to the flux

    [u,g,r,i,z,y]_flux are numpy arrays containing the fluxes
    in ergs/s in each of the LSST bands (0 flux means the star
    is radiating at its quiescent luminosity; 10^32 ergs/s means
    that 10^32 ergs/s need to be added on top of the star's
    quiescent luminosity)
    """

    t_peak_arr, energy_arr = draw_energies(stellar_class, years*365.25, rng)
    duration_arr = duration_from_energy(energy_arr, rng)
    fwhm_arr = fwhm_from_duration(duration_arr)
    del duration_arr
    amplitude_arr = amplitude_from_fwhm_energy(fwhm_arr, energy_arr)
    del energy_arr

    sec_per_year = 365.25*24.0*3600.0
    time_sec_arr = np.arange(0.0, years*sec_per_year, dt)

    johnson_u_flux = np.zeros(len(time_sec_arr))

    fwhm_arr = fwhm_arr*60.0  # convert to seconds
    for amp, fwhm, t_peak in zip(amplitude_arr, fwhm_arr, t_peak_arr*86400.0):

        d_flux = np.piecewise((time_sec_arr-t_peak)/fwhm,
                              [np.logical_and(t_peak-time_sec_arr>=0.0, t_peak-time_sec_arr<=fwhm),
                               np.logical_and(time_sec_arr-t_peak>0.0, time_sec_arr-t_peak<=7.0*fwhm)],
                              [_f_rise_of_t, _f_decay_of_t])

        johnson_u_flux += amp*d_flux

    (u_flux, g_flux, r_flux,
     i_flux, z_flux, y_flux) = lsst_flare_fluxes_from_u(johnson_u_flux)

    return (time_sec_arr/86400.0,
            u_flux, g_flux, r_flux,
            i_flux, z_flux, y_flux)
