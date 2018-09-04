#!/usr/bin/env python3

import numpy as np
import AADG3
from astropy.io import fits
from argparse import ArgumentParser

parser = ArgumentParser(description=
"""Given a FITS output file for a star from S4TESS, this script
rebuilds the AADG3 input files using the data in the headers.  The
namelist, mode data and rotation data will be saved in files with
the basename and extensions .in, .con and .rot""")
parser.add_argument('fits', type=str,
                    help="name of FITS file containing timeseries")
parser.add_argument('basename', type=str,
                    help="basename of output files")
args = parser.parse_args()

f = fits.open(args.fits)
    
nml = {
    'user_seed': f[0].header['SEED'],
    'cadence': 120.0,
    'n_cadences': f[0].header['N_CADS'],
    'n_relax': 4320,
    'n_fine': 50,
    'sig': f[0].header['GRAN_SIG'],
    'rho': 0.45,
    'tau': f[0].header['GRAN_TAU'],
    'inclination': f[0].header['INC'],
    'cycle_period': 100.0,
    'cycle_phase': 0.0,
    'nuac': 0.0,
    'p(1)': 1.52355,
    'p(2)': 0.565349,
    'p(3)': 0.0361707,
    'add_granulation': True,
    'modes_filename': args.basename + '.con',
    'rotation_filename': args.basename + '.rot',
    'output_filename': args.basename + '.asc'
}

AADG3.save_namelist(args.basename + '.in', nml)

modes = f[2].data

with open(args.basename + '.rot', 'w') as f:
    for mode in modes:
        for m in range(1, mode['l']+1):
            f.write('%5i%3i%3i%12.7f\n' % (mode['n'], mode['l'], m, mode['rot']))

modes['rot'] = 0
np.savetxt(args.basename + '.con', modes,
           fmt=['%2i','  %5i','  %12.7e','  %12.7e','  %12.7e','  %12.8e'])

f.close()
