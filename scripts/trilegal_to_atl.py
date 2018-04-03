#!/usr/bin/env python3

import numpy as np
import pandas as pd
from astropy.io import fits
from argparse import ArgumentParser

parser = ArgumentParser(
"""Converts data from TRILEGAL NumPy binary to CSV for ATL."""
)
parser.add_argument('trilegal', type=str,
                    help="name of NumPy binary file containing TRILEGAL data")
parser.add_argument('output', type=str,
                    help="name of output CSV file with input for ATL")
args = parser.parse_args()

# last file is triout_12287.fits

fmts = ['%11.5f','%10.5f','%10.5f','%9.5f','%9.5f','%9.5f','%11.5f']
columns = ('Teff','R','L','gall','galb','imag','d')
header = '%11s %10s %10s %9s %9s %9s %11s\n' % columns

tri = np.load(args.trilegal)
X = np.vstack([10.**tri['logTe'], 10.**tri['logR'],
               10.**tri['logL'], tri['gall'], tri['galb'],
               tri['imag'], 10.**(1.0 + tri['mu0']/5.0)]).T
X = X[tri['comp']>=0]

with open(args.output, 'w') as f:
    f.write(header)
    np.savetxt(f, X, fmt=fmts)

df = pd.DataFrame(X, columns=columns)
df.to_pickle('atl_test.pk', protocol=2)
