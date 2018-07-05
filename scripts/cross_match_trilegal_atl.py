#!/usr/bin/env python

import numpy as np
# import matplotlib.pyplot as pl
from argparse import ArgumentParser

parser = ArgumentParser(description=

"""Cross matches targets in ATL output with original record in
TRILEGAL simulation data, so that we can recover the input for MESA.
Takes as input TRILEGAL data which is reduced to match given ATL data.
Optionally writes a new TRILEGAL-style file that contains the matching
stars."""
                        
)
parser.add_argument('trilegal', type=str,
                    help="name of NumPy binary file containing TRILEGAL data")
parser.add_argument('atl', type=str,
                    help="name of NumPy binary containing ATL output")
parser.add_argument('-o', '--output', type=str, default=None,
                    help="filename for TRILEGAL matches")
parser.add_argument('-v', '--verbose', action='store_true',
                    help="show cross-matching data while script runs")
args = parser.parse_args()

print('Loading TRILEGAL data...')
tri = np.load(args.trilegal)
# print(tri.dtype.names)

print('Loading ATL data...')
atl = np.load(args.atl)
# print(atl.dtype.names)

mu0 = []
chi2 = np.zeros(len(tri))  # []
gall = []
galb = []
chi2 = []
d = 10.**(1.+np.array(tri['mu0'])/5.)

gall = tri['gall']
# gall[gall > 180.] = 360. - gall[gall > 180.]
galb = tri['galb']

print('Cross matching stars...')

if args.verbose:
    print(''.join(['{:>12s}']*6).format('ATL l', 'TRI l',
                                        'ATL b', 'TRI b',
                                        'ATL Teff', 'TRI Teff'))

I = []  # indices of matching stars in TRILEGAL data

for j, row in enumerate(atl):
    allchi2 = (row['GLon']-gall)**2 + (row['GLat']-galb)**2 # + (row['Lum']-10.**tri['logL'])**2/0.01**2
    i = np.argmin(allchi2)
    I.append(i)

    if args.verbose:
        print('{:>12.6f}{:>12.6f}{:>12.6f}{:>12.6f}{:>12.3f}{:>12.3f}'.format(
            row['GLon'], gall[i], row['GLat'], galb[i],
            row['teff'], 10.**tri[i]['logTe']))

if args.output is not None:
    np.save(args.output, tri[I])

