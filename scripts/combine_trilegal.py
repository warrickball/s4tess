#!/usr/bin/env python3

import numpy as np
from astropy.io import fits
from argparse import ArgumentParser

parser = ArgumentParser(description=
                        
"""Reads all the individual TRILEGAL simulation files and combines
them into one NumPy binary file."""

)
parser.add_argument('output', type=str, help="filename for output NumPy binary")
parser.add_argument('--start', type=int, default=0,
                    help="integer index of first TRILEGAL output file (default 0)")
parser.add_argument('--stop', type=int, default=12287,
                    help="integer index of last TRILEGAL output file (default 12287)")
parser.add_argument('--data-dir', type=str,
                    default='/home/wball/rds/miglioa-stellar-grids/TESS_MESA/FitsDir',
                    help="directory containing TRILEGAL output (default=RDS grids)")
parser.add_argument('--basename', type=str, default='triout_{:d}.fits',
                    help="base filename of TRILEGAL output files (default='triout_{:d}.fits')")
parser.add_argument('-v', '--verbose', action='store_const', const=True, default=False,
                    help="print progress")
args = parser.parse_args()

# last file is triout_12287.fits

dtype = [('gall', '>f8'), ('galb', '>f8'), ('Gc', '>i2'), ('logAge', '>f4'), ('M_H', '>f4'), ('m_ini', '>f4'), ('mu0', '>f4'), ('Av', '>f4'), ('comp', '>i2'), ('Mass', '>f4'), ('logL', '>f4'), ('logTe', '>f4'), ('logg', '>f4'), ('label', '>i2'), ('logR', '>f8'), ('logrho', '>f8'), ('numax', '>f4'), ('ddeltanu', '>f4'), ('deltanu', '>f4'), ('pspaci', '>f4'), ('pcoupl', '>f4'), ('mbolmag', '>f4'), ('TESSmag', '>f4'), ('Jmag', '>f4'), ('Hmag', '>f4'), ('Ksmag', '>f4'), ('Keplermag', '>f4'), ('gmag', '>f4'), ('rmag', '>f4'), ('imag', '>f4'), ('zmag', '>f4'), ('DDO51_finfmag', '>f4')]

X = []
for i in range(args.start, args.stop+1):
    filename = args.data_dir + '/' + args.basename.format(i)
    if args.verbose: print('%8i %s' % (i, filename))
    x = fits.getdata(filename, 1, header=False)
    x = x[x['comp'] >= 0]  # individual stars
    x = np.array(x, dtype=dtype)
    X.append(x)

np.save(args.output, np.hstack(X))
