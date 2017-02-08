from __future__ import with_statement
import argparse
import numpy as np
import os
from lsst.sims.catalogs.db import DBObject
from mdwarf_utils import xyz_from_lon_lat_px

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--table", type=str, default=None)
    parser.add_argument("--outdir", type=str, default=None)

    args = parser.parse_args()
    if args.table is None or args.outdir is None:
        raise RuntimeError("must specify table and out dir")

    db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')


    d_z = 50.0
    out_dict = {}

    dtype = np.dtype([('gal_l', float), ('gal_b', float), ('px', float)])

    query = 'SELECT TOP 500000 gal_l, gal_b, parallax FROM %s' % args.table
    chunk_iter = db.get_chunk_iterator(query, chunk_size=100000, dtype=dtype)
    for chunk in chunk_iter:
        xyz = xyz_from_lon_lat_px(chunk['gal_l'], chunk['gal_b'],
                                  0.001*chunk['px'])

        z_quant = np.round(np.abs(xyz[2]/d_z)).astype(int)
        z_dex, z_ct = np.unique(z_quant, return_counts=True)
        for dd, cc in zip(z_dex, z_ct):
            if dd in out_dict:
                out_dict[dd] += cc
            else:
                out_dict[dd] = cc

    with open(os.path.join(args.outdir, "%s_z_dist.txt" % args.table), 'w') as output_file:
        key_list = list(out_dict.keys())
        key_list.sort()
        for kk in key_list:
            output_file.write('%.2f %d\n' % (kk*d_z, out_dict[kk]))

