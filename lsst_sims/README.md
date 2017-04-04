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

#Classifying M-Dwarf stars

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
You can see a direct quotation from his correspondence at the top of
`find_fudge_factor.py`.  Jim recommended scaling the `tau` found for the
"spectroscopically active" population by some fudge factor to get the fraction
of stars that are "flaring active".  The method `_find_fundge_factor` in
`fit_activity_level.py` finds this factor by
fitting the "Active Stars" and "Flaring Stars" curves from Figure 12 of Hilton
et al. 2010 (AJ 140, 1402) to decaying exponentials and taking the ratio of
their respective `tau` parameters.  This script also produces the plot
`plots/fit_fudge_factor.png`, which shows the fits to the two curves.
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
PhD.

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
