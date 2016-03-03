'''
Methods here build the master flare-rate model. This will be used as:

Given the 1 parameter that tunes the flare rate (FFD amplitude),
 what is the cumulative distribution of the relative flux?

I want this relationship be to be analytic, so it will be very fast and
applicable to a zillion sparsely sampled light curves with ease.

'''


import numpy as np
import aflare
import matplotlib.pyplot as plt


def _randomp(num, slope=-2, min=0.1, max=10.):
    # eqn based on this:
    # http://mathworld.wolfram.com/RandomNumber.html
    p = ((max**slope - min**slope) * np.random.random(num) + min**slope)**(1. / slope)
    return p


def SuperLC(ffd_alpha=1, ffd_beta=-1, dur=1):
    '''
    generate a super-sampled (1-minute) light curve of flares for the duration

    FFD fit must be in log(cumulative # / day) _versus_ log(Equiv Dur)

    ffd_beta = slope
    ffd_alpha = intercept
    dur = duration in years
    '''

    dt = 1. / 24. / 60. # 1 minute sampling, in units of days

    time = np.arange(0, dur * 365., dt)
    flux = np.zeros_like(time) # a blank array of fluxes

    # log ED limits:
    ffd_min = -2
    ffd_max = 4

    # to calc the # flares, evaluate the FFD @ the minimum energy,
    # then multiply by the total duration in days
    Nflares = (np.power(10., ffd_min*ffd_beta + ffd_alpha)) * (dur * 365)

    f_energies = _randomp(Nflares, slope=ffd_beta, min=10**ffd_min, max=10**ffd_max)

    # make some REALLY BAD assumptions from event energy to FWHM and Amplitude
    fwhm = (10**((np.log10(f_energies) + 0.5) / 1.5)) / 24. / 60.
    ampl = f_energies/1e1


    # put flares at random places throughout light curve
    t_peak = np.random.random(Nflares) * (dur * 365)

    for k in range(0, int(Nflares)):
        flux = flux + aflare.aflare1(time, t_peak[k], fwhm[k], ampl[k])

    plt.figure()
    plt.plot(time, flux)
    plt.show()

    return flux


def CDist(num=10):
    '''
    For a range of FFD amplitudes (holding the slope fixed)
    - Generate their super-sampled LC's
    - Convert each LC to cumulative fractional flux distribution versus fraction of time
    - Do a 3-D (surface) fit, that parametrizes the cumulative flux distribution versus the 1 free FFD parameter
    - Return the coefficients for this surface fit to be used in the toy model

    '''


    ffd_intercept = np.logspace(-3, 3, num=num)

    for k in range(0, num):
        lc = SuperLC(ffd_b=ffd_intercept[k])
        lc.sort()



    return