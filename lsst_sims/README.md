# Modeling M-Dwarf flare activity

The general strategy we have adopted for simulating M dwarf flaring is to
generate a set of ten-year flaring light curves of different levels of
flaring activity, realistically assign a light curve and a phase offset
to each cool dwarf star in the CatSim simulation, and interpolate the
light curve (incorporating the phase offset) to find the magnitude of
the star at the time of any given simulated observation.

# Classifying levels of flaring activity

In order to correctly assign light curves to the M dwarf stars in our the CatSim
GalFast simulation, we need to divide the stars into five categories:

- Early-type active stars
- Early-type inactive stars
- Mid-type active stars
- Mid-type inactive stars
- Late-type stars (all of which are assumed to be active)

These are the classes into which Eric Hilton's PhD dissertation
(http://search.proquest.com/docview/1069196057) divides
M dwarfs (it is also the classification used by Hilton et al. 2011, 16th
Cambridge Workshop on Cool Stars, Stellar Systems and the Sun, ASP
Conference Series, Vol. 448).  Once we have divided the stars into these
classes, we can use the parameters in Table 4.3 of Eric Hilton's PhD
dissertation to assign flaring light curves to each star.

# Classifying M-Dwarf stars

To find the spectral type of an M-Dwarf star, use the `prob_of_type`
method in `mdwarf_utils.py`.  This method uses the color data from
Table 2 of West et al. 2011 (AJ 141, 97) to calculate the probability
that a given star belongs to one of the classes M0-M9.

Once each star has been assigned a spectral type, we will need to
probabilistically assign each star an activity rate ("active"/"inactive").
The method `find_fraction_spec_active` in the script
`fit_activity_level.py` use the data from
Figure 5 of West et al. 2008 (AJ 135, 785) to model the "active"
fraction of each spectral type as a function of distance from the
Galactic plane with a decaying exponential of the form `A*exp[-z/tau] + B`
where `z` is the distance in parsecs from the Galactic plane and `A`, `tau`,
and `B` are parameters to be fit.  Running `fit_activity_level.py` directly
will produce the plot `plots/exp_decay_by_type.png` showing the
fits to the West et al. data.

The fitting done by `find_fraction_spec_active` actually returns the fraction
of each spectral type that is "spectroscopically active."
According to Jim Davenport,
there is a difference between "spectroscopically active" and "flaring active".
Quoting a private correspondence from him:

    If I remember correctly what I said, it's that spectroscopic
    activity doesn't == flaring activity. So stars could be "magnetically active"
    as determined by Halpha emission in their spectra, but not "flare active"
    (not as well defined)

    The evidence of this is Hilton et al (2010),
    http://adsabs.harvard.edu/abs/2010AJ....140.1402H (Fig 12) where they
    find the scale height of flare stars is lower than for "active" stars.

    Of course as the specialist I have to say "this is probably more
    complicated"...

    So my off-the-cuff advice:
    do exactly as you propose, use West et al. to make "active vs inactive"
    as a function of height, assign flares, etc. BUT, adopt some correction
    (or ignorance) scaling factor that makes the scale height lower. I don't
    know what that # should be, but you could extract something plausible
    out of Fig 12 from Hilton2010 (e.g. 80pc instead of 120pc)

Jim recommended scaling the `tau` found for the
"spectroscopically active" population by some fudge factor to get the fraction
of stars that are "flaring active".  The method `_find_fudge_factor` in
`fit_activity_level.py` finds this factor by
fitting the "Active Stars" and "Flaring Stars" curves from Figure 12 of Hilton
et al. 2010 (AJ 140, 1402) to decaying exponentials and taking the ratio of
their respective `tau` parameters.
The method `find_fraction_flare_active` uses `find_fraction_spec_active` and
`_find_fudge_factor` to determine the fraction of each spectral type (M0-M8)
that is flaring active as a function of distance from the Galactic Plane.

The method `activity_type_from_color_z` reads in a star (or stars') r-i and i-z
colors and its distance from the Galactic Plane and, using a random number
generator from numpy as well as all of the infrastructure described in above,
assigns the star to one of the variability classes mentioned at the top of this
document (early active, early inactive, mid active, mid inactive, late active).

# Generating light curves

We have elected to pre-generate four ten-year light curves for each of the
flaring activity classes ('early active', 'early inactive', 'mid active', etc.)
in our simulation.  The method `light_curve_from_class` in `mdwarf_utils.py`
takes an activity class, a number of years, and a numpy random number
generator and returns a light curve in all six LSST bands simulated over
the specified number of years.  Below, we describe the steps taken by
`light_curve_from_class`.

To determine the number and energy of flares experienced by a star, we use the
cumulative flare distributions as a function of energy presented in Table 4.3 of
Eric Hilton's PhD dissertation (http://search.proquest.com/docview/1069196057).
The method `draw_energies` in `mdwarf_utils.py` takes a flaring activity class,
a duration in days, and a numpy random number generator, and randomly draws
flares according to distributions specified by Hilton.  This method will return
the peak times of the flares in days (randomly drawn from a uniform distribution
spread out over the specified duration) and the total energy of the flares in
ergs in the Johnson U band (randomly drawn from the Hilton distributions with a
maximum flare energy set at 10^34 ergs).  The script `validate_energy.py`
validates this process by drawing 10 years worth of flares and plotting the
simulated cumulative flare rate against the actual distributions from Hilton's
PhD.  This plot will be found in `plots/energy_dist_validation.png` after
running `validate_energy.py`.

Once we have a list of flares, their time of peak, and their total energies, we
need to model their actual rise and fall.  We do this using the profile
presented in Davenport et al. 2014 (ApJ 797, 122) equations (1) and (4).  This
profiles requires as input the time full width half maximum of the flare (i.e.
the width of the profile in time between rising to one half of its maximum flux
and falling back to one half of its maximum flux) and an amplitude (i.e. the
maxmimum flux of the flare).  We model the relationship between energy, FWHM,
and amplitude using observations of the well-known flaring star GJ 1243
presented in Hawley et al. 2014 (ApJ 797, 121) and provided by Jim Davenport at
(http://github.com/jradavenport/GJ1243-Flares).

To find the FWHM time duration of the flares, we bin the duration versus energy
plot from Hawley et al. 2014 (middle panel of Figure 10) into energy bins of 1.0
dex, and fit the duration distribution within each energy bin as a Gaussian.  To
find the durations of our simulated flares, we draw randomly from these
Gaussians, based on the energy of the flare.  This is done in the method
`duration_from_energy` in `mdwarf_utils.py`.  The script `validate_duration.py`
validates this process by taking the actual energies of the Hawley et al. 2014
flares, feeding them through `duration_from_energy`, and plotting the resulting
distribution of durations against the actual durations in the data.  These plots
will be found in `plots/duration_plot.png` after running `validate_duration.py`.

Once we have flare energies and FWHM time durations, the amplitude is fixed. The
method `amplitude_from_fwhm_energy` in `mdwarf_utils.py` takes a list of FWHM
time durations and energies in the Johnson U band and returns the amplitude of
the flare required by the Davenport et al. 2014 profile and the energy.  The
script `validate_amplitude.py` validates this process by taking the actual
energies and durations reported in the Hawley et al. 2014 data and running them
through `amplitude_from_fwhm_energy` and plotting the resulting distribution of
amplitudes versus the actual amplitudes reported in the data.  These plots will
be found in `plots/amplitude_plot.png` after running `validate_amplitude.py`.

The entire process as described above is validated by `validate_end_to_end.py`.
This script simulates flares on a `mid active` star (like GJ 1243) covering a
duration of time equal to the total duration of time covered by the Hawley et
al. 2014 data.  These simulated flares are then run through
`duration_from_energy` and `amplitude_from_duration_energy`.  The resulting
distributions of amplitude, energy, and duration are plotted against the actual
datain `plots/end_to_end_1D.png` and `plots/end_to_end_2D.png`.

At this point, `light_curve_from_class` has all of the information it needs to
create a flaring light curve from the star in the Johnson U band.  To convert to
LSST bands, we model all flares as 9000K black bodies.  The method
`lsst_flare_fluxes_from_u` in `mdwarf_utils.py` finds the scaling factor
by which we must multiply a fiducial blackbody curve to get its total U band
flux to equal the U band flux of the flare.  `lsst_flare_fluxes_from_u` then
applies that same factor to the LSST band fluxes of the fiducial black body
curve.  This gives us flaring light curve fluxes in all six LSST bands.  Note:
these light curves represent the total energy output by the star in all
directions over the sky.  In order to convert to observed fluxes, you must
multiply by (effective area)/(4*pi*R^2) where R is the distance to the star.
The script `validate_magnitudes.py` queries one star from fatboy in each flaring
variability class and outputs the flaring delta magnitudes into the file
`delta_m_data.txt`.  This data can be compared to Figure 10 of Chang et al. 2015
(ApJ 814, 35) for validation purposes.
