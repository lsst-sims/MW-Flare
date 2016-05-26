import numpy as np
import matplotlib.pyplot as plt


def generate_visits(Nvisits=900, tspan=10, stat=False,
                    seasonscale=365./5):
    '''
    Use some very crude approximations for how visits will be spaced out:

    - Survey starts at midnight, time = 0.0
    - Can only observe at night, time > 0.75 | time < 0.25
    - Exposures are clustered around a season w/ a gaussian shape each year
    - Field is observable for first half of year, 0 < date < 182
    - On average, each field should be hit every 3 days during observable season

    Set "stat=True" if you want a plot and a couple basic statistics about the cadence

    '''

    # generate random times for visit, between [0.75 and 0.25]
    time_of_day = np.random.random(Nvisits)/2. - 0.25

    date_of_year = np.floor(np.random.normal(loc=365./4., scale=seasonscale, size=Nvisits))

    year_of_obs = np.floor(np.random.random(Nvisits) * tspan) * 365.

    date_obs = time_of_day + date_of_year + year_of_obs

    date_obs.sort()

    if stat is True:
        print('mean time between visits:')
        print(np.mean(date_obs[1:] - date_obs[:-1]))

        print('median time between visits:')
        print(np.median(date_obs[1:] - date_obs[:-1]))

        plt.figure()
        _ = plt.hist(date_obs, bins=np.arange(date_obs.min(), date_obs.max(),7),
                     histtype='stepfilled', color='k')
        plt.xlabel('Time (days)')
        plt.ylabel('# Visits per Week')
        plt.show()

    return date_obs


def photerror(mag, nfloor=0.005, m5=24.89, gamma=0.038):
    '''
    http://arxiv.org/pdf/0805.2366v4.pdf
    use values from Ivezic 2008, assuming g-band

    gamma is very insensitive to filter
    '''

    x = 10.0**(0.4*(mag - m5))
    sigrand2 = (0.04 - gamma) * x + gamma * (x**2.0)

    err = (nfloor**2. + sigrand2)**0.5
    return err


def downsample(time,flux):
    '''
    take super-sampled LC (from flare_prob), uses simple linear interpretation
    to down-sample to LSST cadence.

    Assumes 10 years of LSST with 900 visits
    '''

    tout = generate_visits()

    fout = np.interp(tout, time, flux)

    return tout, fout


def opsim_cadence(time, flux, mag, los, random=True):
    '''
    Instead of my LSST Toy Model for cadence, use the OpSim inputs Ruth provided.
    Assume u- and g-band are usable

    los = the latitude, from -10 to -80
    time, flux = the super-sampled Kepler-like light curve
    mag = the apparent mag of the target star

    random = True means pick random places within the super-sampled light curve
        Use this because I only made 6 super-sampled light curves, and every star will have
        same flare properties based on timestamp
    '''

    dir = 'cadence_files/'
    ufile = 'l45b' + str(int(los)) + '.0_u_cadence.txt'
    gfile = 'l45b' + str(int(los)) + '.0_g_cadence.txt'

    tu, mu = np.loadtxt(dir + ufile, dtype='float', unpack=True)
    tg, mg = np.loadtxt(dir + gfile, dtype='float', unpack=True)

    u2g = 2.0  # ratio of u-band to g-band flare flux (ad hoc)
    sclu = np.ones_like(tu) * u2g
    sclg = np.ones_like(tg)

    tout = np.concatenate((tu, tg))
    m5 = np.concatenate((mu, mg))
    scl = np.concatenate((sclu, sclg))

    ss = np.argsort(tout)
    tout = tout[ss]
    m5 = m5[ss]
    scl = scl[ss]

    if random is True:
        time2 = np.random.choice(time, len(time))
    else:
        time2=time

    fout = np.interp(tout, time2, flux)

    err = photerror(mag, m5=m5)

    return tout, fout, err, scl


if __name__ == "__main__":
    generate_visits(stat=True)
