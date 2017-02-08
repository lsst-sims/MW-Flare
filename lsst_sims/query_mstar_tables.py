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

    parser.add_argument('--max_type', type=int, help='maximum M sub-type',
                        default=8)

    args = parser.parse_args()
    if args.suffix is None:
        raise RuntimeError("Must specify a suffix")

    # from Figure 2 of Kowalski et al 2009 (AJ 138, 633)
    # anything outside these bounds is considered later than M-dwarf
    r_i_cutoff = 2.49 + 3*0.027
    i_z_cutoff = 1.35 + 3*0.0095
    d_color = 0.01
    n_later = 0
    n_total = 0

    # meant to correspond with figure 4 of West et al. 2008
    # (AJ 135, 785)
    z_bins = [(ii-25.0, ii) for ii in np.arange(25.0, 250.0, 25.0)]

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    color_color_grid = {}
    star_counts_zbin = {}
    for bin in z_bins:
        star_counts_zbin[bin] = {}
        for ix in range(args.max_type+1):
            star_counts_zbin[bin]['M%d' % ix] = 0
        star_counts_zbin[bin]['later'] = 0

    table_name = 'stars_mlt_part_%s' % args.suffix
    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')

    dtype = np.dtype([('lon', float), ('lat', float), ('px', float),
                      ('u', float), ('g', float), ('r', float),
                      ('i', float), ('z', float)])

    query = 'SELECT gal_l, gal_b, parallax, sdssu, sdssg, sdssr, sdssi, sdssz '
    query += 'FROM %s' % table_name

    chunk_iter = db.get_chunk_iterator(query, chunk_size=100000, dtype=dtype)

    for m_stars in chunk_iter:

        n_total += len(m_stars)

        ri = m_stars['r'] - m_stars['i']
        iz = m_stars['i'] - m_stars['z']
        ri_dex = (ri/d_color).astype(int)
        iz_dex = (iz/d_color).astype(int)
        assert iz_dex.max() < 1000
        color_color_dex = ri_dex*10000+iz_dex
        color_color_unique, color_color_counts = np.unique(color_color_dex,
                                                           return_counts=True)

        ri_fin = color_color_unique//10000
        iz_fin = color_color_unique%10000

        for ri, iz, ct in zip(ri_fin, iz_fin, color_color_counts):
            if ri in color_color_grid:
                if iz in color_color_grid[ri]:
                    color_color_grid[ri][iz] += ct
                else:
                    color_color_grid[ri][iz] = ct
            else:
                color_color_grid[ri] = {}
                color_color_grid[ri][iz] = ct

        xyz = xyz_from_lon_lat_px(m_stars['lon'], m_stars['lat'],
                                  0.001*m_stars['px'])

        for bin in z_bins:
            local_dexes = np.where(np.logical_and(np.abs(xyz[2])>bin[0],
                                                  np.abs(xyz[2])<=bin[1]))

            local_m_stars = m_stars[local_dexes]
            good_colors = np.where(np.logical_and(local_m_stars['r']-local_m_stars['i']<r_i_cutoff,
                                                  local_m_stars['i']-local_m_stars['z']<i_z_cutoff))

            local_later = len(local_m_stars) - len(good_colors[0])
            n_later += local_later
            actual_m_stars = local_m_stars[good_colors]

            prob = prob_of_type(actual_m_stars['r']-actual_m_stars['i'],
                                actual_m_stars['i']-actual_m_stars['z']).transpose()

            assert prob.shape[0] == len(actual_m_stars)

            local_types = np.argmax(prob, axis=1)
            assert len(local_types) == len(actual_m_stars)

            unique_types, unique_counts = np.unique(local_types,
                                                    return_counts=True)
            for tt, cc in zip(unique_types, unique_counts):
                star_counts_zbin[bin]['M%d' % tt] += cc
            star_counts_zbin[bin]['later'] += local_later

    for bin in z_bins:
        out_name = os.path.join(args.out_dir,
                                'mdwarf_count_%s_%d_%d.txt' %
                                (args.suffix, bin[0], bin[1]))

        with open(out_name, 'w') as output_file:
            for ix in range(args.max_type+1):
                output_file.write('M%d: %d\n' % (ix, star_counts_zbin[bin]['M%d' % ix]))
            output_file.write('later: %d\n' % star_counts_zbin[bin]['later'])

    out_name = os.path.join(args.out_dir,
                            'color_color_grid_%s.txt' % args.suffix)
    with open(out_name, 'w') as output_file:
        output_file.write('# r-i i-z ct\n')
        ri_list = list(color_color_grid.keys())
        ri_list.sort()
        for ri in ri_list:
            iz_list = list(color_color_grid[ri].keys())
            iz_list.sort()
            for iz in iz_list:
                output_file.write('%.2f %.2f %d\n' % (ri*d_color, iz*d_color,
                                                      color_color_grid[ri][iz]))

    print("n_later %d of %d\n" % (n_later, n_total))
