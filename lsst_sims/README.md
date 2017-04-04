#Modeling M-Dwarf flare activity

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
