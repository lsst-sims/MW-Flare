from __future__ import with_statement
import os
from lsst.utils import getPackageDir
from lsst.sims.photUtils import BandpassDict, Sed

# load the baseline LSST bandpasses
bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

mlt_dir = os.path.join(getPackageDir('sims_sed_library'),
                       'starSED', 'mlt')

mlt_list = os.listdir(mlt_dir)

with open('mlt_spectra_mags.txt', 'w') as output:
    output.write('# name u g r i z y\n')
    for mlt_name in mlt_list:
        full_name = os.path.join(mlt_dir, mlt_name)
        ss = Sed()
        ss.readSED_flambda(full_name)
        mag_dict = bp_dict.magDictForSed(ss)
        output.write('%s %.6g %.6g %.6g %.6g %.6g %.6g\n'
                     % (mlt_name, mag_dict['u'], mag_dict['g'],
                        mag_dict['r'], mag_dict['i'], mag_dict['z'],
                        mag_dict['y']))
