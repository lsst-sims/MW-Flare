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
The script `fit_by_type.py` uses the data from Figure 5 of West et al. 2008
(AJ 135, 785) to model the "active" fraction of each spectral type as a
function of its distances from the Galactic plane with a decaying exponential
of the form `A*exp[-z/tau] + B` where `z` is the distance in parsecs from the
Galactic plane and `A`, `tau`, and `B` are parameters to be fit.  The resulting
fits are saved to text files in the `type_fits/` directory of this repository.
This script will also produce the plot `plots/exp_decay_by_type.png` showing the
fits to the West et al. data.

If you would like to validate the results of this fit, run the script
`query_mstar_tables.py` on each of the MLT dwarf star partitions of the CatSim
database (0870, 1100, 1160, 1180, 1200, 1220, 1250, 1400), making sure the
output is sent to the `query_results/` directory.  Then run
`plot_cumulative_activity.py`, which will compare the cumulative fraction of
active stars found by our model with that found in Hilton et al. 2010 (AJ 140,
1402).  It will also compare the total fraction of each spectral class that is
active with West et al Figure 3.

The fitting done by `fit_by_type.py` actually returns the fraction of each
spectral type that is "spectroscopically active."  According to Jim Davenport,
there is a difference between "spectroscopically active" and "flaring active".
You can see a direct quotation from his correspondence at the top of
`find_fudge_factor.py`.  Jim recommended scaling the `tau` found for the
"spectroscopically active" population by some fudge factor to get the fraction
of stars that are "flaring active".  `find_fudge_factor.py` finds this factor by
fitting the "Active Stars" and "Flaring Stars" curves from Figure 12 of Hilton
et al. 2010 (AJ 140, 1402) to decaying exponentials and taking the ratio of
their respective `tau` parameters.  This script also produces the plot
`plots/fit_fudge_factor.png`, which shows the fits to the two curves.

We should be able to assign each star a probability of being flaring
active/inactive by finding its spectral type, finding its distance above the
Galactic plane (using the method `xyz_from_lon_lat_px` defined in
`mdwarf_utils.py`), and using the fit from `fit_by_type.py` (adjusting `tau` by
the fudge factor found above).
