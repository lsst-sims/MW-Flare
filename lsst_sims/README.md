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
