from __future__ import with_statement
import argparse
import os
import sys
import numpy as np
from lsst.sims.catalogs.db import DBObject
from mdwarf_utils import prob_of_type, xyz_from_lon_lat_px

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='get vital stats from an '
                                                'mdwarf table on fatboy')

    parser.add_argument('--suffix', type=str, help='suffix of table name',
                        default=None)

    parser.add_argument('--out_dir', type=str, help='output directory',
                        default='output')

    args = parser.parse_args()
    if args.suffix is None:
        raise RuntimeError("Must specify a suffix")

    # from Figure 2 of Kowalski et al 2009 (AJ 138, 633)
    # anything outside these bounds is considered later than M-dwarf
    r_i_cutoff = 2.3
    i_z_cutoff = 1.35
    n_later = 0
    n_total = 0

    # meant to correspond with figure 4 of West et al. 2008
    # (AJ 135, 785)
    z_bins = [(ii-25.0, ii) for ii in np.arange(25.0, 250.0, 25.0)]

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    star_counts = {}
    for bin in z_bins:
        star_counts[bin] = {}
        for ix in range(8):
            star_counts[bin]['M%d' % ix] = 0

    table_name = 'stars_mlt_part_%s' % args.suffix
    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

    dtype = np.dtype([('lon', float), ('lat', float), ('px', float),
                      ('u', float), ('g', float), ('r', float),
                      ('i', float), ('z', float)])

    query = 'SELECT gal_l, gal_b, parallax, sdssu, sdssg, sdssr, sdssi, sdssz '
    query += 'FROM %s' % table_name

    chunk_iter = db.get_chunk_iterator(query, chunk_size=100000, dtype=dtype)

    for chunk in chunk_iter:
        good_dexes = np.where(np.logical_and(chunk['r']-chunk['i']<r_i_cutoff,
                                             chunk['i']-chunk['z']<i_z_cutoff))


        n_later += len(chunk)-len(good_dexes[0])
        n_total += len(chunk)
        print('n_later: %.3e of %.3e' % (n_later, n_total))

        m_stars = chunk[good_dexes]
        xyz = xyz_from_lon_lat_px(m_stars['lon'], m_stars['lat'],
                                  0.001*m_stars['px'])

        for bin in z_bins:
            local_dexes = np.where(np.logical_and(np.abs(xyz[2])>bin[0],
                                                  np.abs(xyz[2])<=bin[1]))

            local_m_stars = m_stars[local_dexes]
            prob = prob_of_type(local_m_stars['r']-local_m_stars['i'],
                                local_m_stars['i']-local_m_stars['z']).transpose()

            assert prob.shape[0] == len(local_m_stars)
            local_types = np.argmax(prob, axis=1)
            assert len(local_types) == len(local_m_stars)

            unique_types, unique_counts = np.unique(local_types,
                                                    return_counts=True)
            for tt, cc in zip(unique_types, unique_counts):
                star_counts[bin]['M%d' % tt] += cc

    for bin in z_bins:
        out_name = os.path.join(args.out_dir,
                                'mdwarf_count_%s_%d_%d.txt' %
                                (args.suffix, bin[0], bin[1]))

        with open(out_name, 'w') as output_file:
            for ix in range(8):
                output_file.write('M%d: %d\n' % (ix, star_counts[bin]['M%d' % ix]))

    print("n_later %d of %d\n" % (n_later, n_total))
