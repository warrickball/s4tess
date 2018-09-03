#!/usr/bin/env python

import pandas as pd
from astropy.io import fits
from argparse import ArgumentParser

parser = ArgumentParser(description=
"""Takes a list of FITS files and writes their headers to a CSV file.
By default, removes data for FITS defaults 'SIMPLE', 'BITPIX',
'NAXIS', 'EXTEND', 'ORIGIN' and 'DATE', if present.""")
parser.add_argument('table', type=str,
                    help="name for output CSV file")
parser.add_argument('fits_filenames', type=str, nargs='+',
                    help="name of FITS files from which to construct "
                    "table, in the desired order")
parser.add_argument('--remove', type=str, nargs='+',
                    default=['SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'ORIGIN', 'DATE'],
                    help="list of header keys to remove "
                    "(default=['SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'ORIGIN', 'DATE'])")
args = parser.parse_args()

# read one file to get the columns in order
columns = [k for k, v in fits.getheader(args.fits_filenames[0]).items()]
for k in args.remove:
    try:
        columns.remove(k)
    except ValueError:
        continue
        
pd.DataFrame([{k:v for k,v in fits.getheader(filename).items()}
              for filename in args.fits_filenames]).to_csv(args.table, columns=columns, index=False)
