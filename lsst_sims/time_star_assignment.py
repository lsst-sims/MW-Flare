from mdwarf_utils import xyz_from_lon_lat_px
from fit_activity_level import find_fraction_flare_active

import numpy as np

import time

n_stars = 100000
rng = np.random.RandomState(66)
lon = rng.random_sample(n_stars)*360.0
lat = rng.random_sample(n_stars)*180.0 - 90.0
parallax = rng.random_sample(n_stars)*0.1

class_poss = ['M%d' % ii for ii in range(9)]

class_poss = np.array(class_poss)

class_list = class_poss[rng.randint(0,9,n_stars)]

t_start = time.time()
ff_control = []
for ll, bb, px, st in zip(lon, lat, parallax, class_list):
    xyz = xyz_from_lon_lat_px(ll, bb, px)
    ff = find_fraction_flare_active(st, xyz[2])
    ff_control.append(ff)
print '%d stars took %e\n' % (n_stars, time.time()-t_start)

t_start = time.time()
xyz = xyz_from_lon_lat_px(lon, lat, parallax)
ff_test = find_fraction_flare_active(class_list, xyz[2])
print 'took %e vectorized' % (time.time()-t_start)

np.testing.assert_array_equal(ff_control, ff_test)

r_i = rng.random_sample(n_stars)*3.0
i_z = rng.random_sample(n_stars)*2.0

from mdwarf_utils import activity_type_from_color_z

t_start = time.time()
xyz = xyz_from_lon_lat_px(lon, lat, parallax)
type_list = activity_type_from_color_z(r_i, i_z, xyz[2], rng)
print 'whole process took %e on %d' % (time.time()-t_start,n_stars)

