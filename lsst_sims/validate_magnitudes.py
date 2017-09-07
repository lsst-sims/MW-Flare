"""
This script will produce delta_magnitude scatter plots of simulated flares
to compare with Figures 10 and 14 of Chang et al. 2015 (ApJ 814, 35)
"""
from __future__ import with_statement
from lsst.sims.catalogs.db import DBObject
import numpy as np

rng = np.random.RandomState(813)
_au_to_parsec = 1.0/206265.0

table = 'stars_mlt_part_1180'

db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
              port=1433, driver='mssql+pymssql')

query = 'SELECT TOP 100 gal_l, gal_b, parallax, '
query += 'sdssr, sdssi, sdssz, '
query += 'umag, gmag, rmag, imag, zmag, ymag '
query += 'FROM %s' % table

dtype = np.dtype([('lon', float), ('lat', float), ('parallax', float),
                  ('sdssr', float), ('sdssi', float), ('sdssz', float),
                  ('umag', float), ('gmag', float), ('rmag', float),
                  ('imag', float), ('zmag', float), ('ymag', float)])

data = db.execute_arbitrary(query, dtype=dtype)

# parallax will be in milli arcsec
# magnorm = -2.5*log10(flux_scale)-18.402732642

from mdwarf_utils import xyz_from_lon_lat_px

xyz = xyz_from_lon_lat_px(np.degrees(data['lon']), np.degrees(data['lat']),
                          0.001*data['parallax'])

from mdwarf_utils import prob_of_type
prob = prob_of_type(data['sdssr']-data['sdssi'], data['sdssi']-data['sdssz'])

import time
t_start = time.time()
prob_i_vec = [np.argmax(p_row) for p_row in prob.transpose()]
prob_vec = ['M%d' % ii if ii<9 else 'M8' for ii in prob_i_vec]

from fit_activity_level import find_fraction_flare_active

frac = find_fraction_flare_active(prob_vec, np.abs(xyz[2]))
draws = rng.random_sample(len(frac))
activity_level = ['active' if dd<ff else 'inactive'
                  for ff, dd in zip(frac, draws)]

type_dict = {}
for ii in range(3):
    type_dict['M%d' % ii] = 'early'
for ii in range(3,6):
    type_dict['M%d' % ii] = 'mid'
for ii in range(6,9):
    type_dict['M%d' % ii] = 'late'

from lsst.sims.utils import radiansFromArcsec
dd = _au_to_parsec/radiansFromArcsec(data['parallax']*0.001)
from lsst.sims.photUtils import PhotometricParameters

params = PhotometricParameters()

# divide the effective area of the LSST mirror by the sphere
# through which the flare is radiating (3.08576e18 is a parsec
# in cm)
flux_factor = params.effarea/(4.0*np.pi*np.power(dd*3.08576e18, 2))

types_simulated = []
du = []
dg = []
dr = []
di = []
dz = []
dy = []

from mdwarf_utils import light_curve_from_class
from lsst.sims.photUtils import Sed

ss = Sed()
base_u = ss.fluxFromMag(data['umag'])
base_g = ss.fluxFromMag(data['gmag'])
base_r = ss.fluxFromMag(data['rmag'])
base_i = ss.fluxFromMag(data['imag'])
base_z = ss.fluxFromMag(data['zmag'])
base_y = ss.fluxFromMag(data['ymag'])

for i_star in range(len(data)):
    thermo_class = type_dict[prob_vec[i_star]]
    if thermo_class == 'late':
        activity = 'active'
    else:
        activity = activity_level[i_star]

    flare_type = '%s_%s' % (thermo_class, activity)
    if flare_type in types_simulated:
        continue
    types_simulated.append(flare_type)
    print 'simulating ',flare_type
    (time,
     uflare, gflare, rflare,
     iflare, zflare, yflare,
     tpeak) = light_curve_from_class(flare_type, 2.0, rng)


    uflare *= flux_factor[i_star]
    gflare *= flux_factor[i_star]
    rflare *= flux_factor[i_star]
    iflare *= flux_factor[i_star]
    zflare *= flux_factor[i_star]
    yflare *= flux_factor[i_star]

    uflux = np.interp(tpeak, time, uflare)
    gflux = np.interp(tpeak, time, gflare)
    rflux = np.interp(tpeak, time, rflare)
    iflux = np.interp(tpeak, time, iflare)
    zflux = np.interp(tpeak, time, zflare)
    yflux = np.interp(tpeak, time, yflare)

    du += list(data['umag'][i_star] - ss.magFromFlux(base_u[i_star]+uflux))
    dg += list(data['gmag'][i_star] - ss.magFromFlux(base_g[i_star]+gflux))
    dr += list(data['rmag'][i_star] - ss.magFromFlux(base_r[i_star]+rflux))
    di += list(data['imag'][i_star] - ss.magFromFlux(base_i[i_star]+iflux))
    dz += list(data['zmag'][i_star] - ss.magFromFlux(base_z[i_star]+zflux))
    dy += list(data['ymag'][i_star] - ss.magFromFlux(base_y[i_star]+yflux))


with open('delta_m_data.txt', 'w') as output_file:
    for uu, gg, rr, ii, zz, yy in zip(du, dg, dr, di, dz, dy):
        output_file.write('%e %e %e %e %e %e\n' % (uu, gg, rr, ii, zz, yy))
