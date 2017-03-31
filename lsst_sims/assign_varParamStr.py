from __future__ import print_function
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--table', type=str, default=None,
                    help='The name of the table on fatboy that you want '
                          'to process')
parser.add_argument('--out_dir', type=str, default=None,
                    help='The directory where we will write the output of '
                         'this script')
parser.add_argument('--use_tunnel', type=bool, default=False,
                    help='Do we need to use the SSH tunnel (because we '
                         'are running from a machine that is not on '
                         'UW campus)?')
parser.add_argument('--limit', type=int, default=None,
                    help='Maximum number of stars to process '
                         '(for testing)')
parser.add_argument('--seed', type=int, default=None,
                     help='Seed for random number generator')
parser.add_argument('--n_curves', type=int, default=4,
                    help='The number of flares to draw from for each '
                         'level of flaring activity')
parser.add_argument('--chunk_size', type=int, default=100000,
                    help='How many stars to download from the '
                         'database at once')

args = parser.parse_args()

if args.table is None:
    raise RuntimeError("must specify a table")
if args.out_dir is None:
    raise RuntimeError("must specify out_dir")
if args.seed is None:
    raise RuntimeError("must specify a seed")

import os

if not os.path.isdir(args.out_dir) and not os.path.exists(args.out_dir):
    os.mkdir(args.out_dir)
if not os.path.isdir(args.out_dir):
    raise RuntimeError("%s is not a directory" % args.out_dir)

from lsst.sims.catalogs.db import DBObject

if args.use_tunnel:
    from lsst.sims.catUtils.baseCatalogModel import BaseCatalogConfig
    config = BaseCatalogConfig
    host = config.host
    port = config.port
    database = config.database
    driver = config.driver
else:
    host = 'fatboy.phys.washington.edu'
    port = 1433
    database = 'LSSTCATSIM'
    driver = 'mssql+pymssql'

db = DBObject(database=database, host=host, port=port, driver=driver)

if args.limit is None:
    query = 'SELECT '
else:
    query = 'SELECT TOP %d ' % args.limit

query += 'simobjid, gal_l, gal_b, parallax, sdssr, sdssi, sdssz '
query += 'FROM %s ' % args.table

from mdwarf_utils import activity_type_from_color_z
from mdwarf_utils import xyz_from_lon_lat_px
import os
import json
import numpy as np
import time

t_start = time.time()

rng = np.random.RandomState(args.seed)

dtype = np.dtype([('id', int),
                  ('lon', float), ('lat', float), ('parallax', float),
                  ('sdssr', float), ('sdssi', float), ('sdssz', float)])

chunk_iter = db.get_arbitrary_chunk_iterator(query, chunk_size=args.chunk_size,
                                             dtype=dtype)

out_name = os.path.join(args.out_dir, '%s_flaring_varParamStr.txt' % args.table)
print('outname ',out_name)
has_written = False

for data_chunk in chunk_iter:
    xyz = xyz_from_lon_lat_px(np.degrees(data_chunk['lon']),
                              np.degrees(data_chunk['lat']),
                              data_chunk['parallax']*0.001)

    (activity_class,
     spectral_type) = activity_type_from_color_z(data_chunk['sdssr']-data_chunk['sdssi'],
                                                 data_chunk['sdssi']-data_chunk['sdssz'],
                                                 xyz[2], rng)

    lc_indices = rng.randint(0, args.n_curves, len(data_chunk))
    offset = rng.random_sample(len(data_chunk))*3652.5
    varParamStr_list = [json.dumps({'pars':{'lc':'%s_%d.txt' % (aa, ii), 't0': oo},
                                    'varMethodName': 'applyMLTflaring'})
                        for aa,ii,oo in zip(activity_class, lc_indices,offset)]

    len_str = [len(vv) for vv in varParamStr_list]
    print('max len ',max(len_str))

    if has_written:
        status = 'a'
    else:
        status = 'w'
        has_written = True

    with open(out_name, status) as out_file:
        for ii, vv, mm, zz in zip(data_chunk['id'], varParamStr_list,
                                  spectral_type, xyz[2]):
            out_file.write("%d; %s; %s; %e\n" % (ii, vv, mm, zz))

print('that took %e seconds' % (time.time()-t_start))
