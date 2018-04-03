#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from argparse import ArgumentParser

parser = ArgumentParser(description=

"""Cross matches targets in ATL output with original record in
TRILEGAL simulation data, so that we can recover the input for
MESA."""
                        
)
parser.add_argument('--trilegal', '--TRILEGAL', type=str, default='../data/tri.npy',
                    help="name of NumPy binary file containing TRILEGAL data")
parser.add_argument('--atl', '--ATL', type=str, default='../data/atl.npy',
                    help="name of NumPy binary containing ATL output")
parser.add_argument('--csv-fallback', type=str,
                    default='/home/wball/mnt/adfbison/data/trilegal_results.csv',
                    help="name of fallback CSV file containing ATL output")
args = parser.parse_args()

print('Loading TRILEGAL data...')
tri = np.load(args.trilegal)
    
print(tri.dtype.names)

print('Loading ATL data...')
try:
    atl = np.load(args.atl)
except IOError:
    print("Falling back to original CSV...")
    atl = np.genfromtxt(args.csv_fallback, delimiter=',', names=True)
    np.save(args.atl, atl)
    
print(atl.dtype.names)

mu0 = []
chi2 = np.zeros(len(tri))  # []
gall = []
galb = []
chi2 = []
I = np.zeros(len(atl), dtype=int) # []
d = 10.**(1.+np.array(tri['mu0'])/5.)

gall = tri['gall']
gall[gall > 180.] = 360. - gall[gall > 180.]
galb = tri['galb']

print(''.join(['{:>12s}']*6).format('ATL l', 'TRI l', 'ATL b', 'TRI b',
                                    'ATL Teff', 'TRI Teff'))

# for j, (GLon, GLat) in enumerate(zip(atl[:10]['GLon'], atl[:10]['GLat'])):
for j, row in enumerate(atl[:10]):
    allchi2 = (row['GLon']-gall)**2 + (row['GLat']-galb)**2 # + (row['Lum']-10.**tri['logL'])**2/0.01**2
    # allchi2 = (GLon-gall)**2 + (GLat-galb)**2 # + (row['Lum']-10.**tri['logL'])**2/0.01**2
    # allchi2 = (row['d']-d)**2/0.01**2 + (row['teff']-10.**tri['logTe'])**2 + (row['Lum']-10.**tri['logL'])**2/0.01**2
    i = np.argmin(allchi2)
    # chi2.append(allchi2[i])
    I[j] = i  # I.append(i)
    
    # if tri[i]['gall'] < 180.0:
    #     gall = tri[i]['gall']
    # else:
    #     gall = 360.0 - tri[i]['gall']
        
    print('{:>12.6f}{:>12.6f}{:>12.6f}{:>12.6f}{:>12.3f}{:>12.3f}'.format(
        row['GLon'], gall[i], row['GLat'], galb[i],
        row['teff'], 10.**tri[i]['logTe']))

np.save('match_index.npy', I)

# pl.semilogy(tri['gall'][I], chi2, 'o')
# pl.xlabel('gal. longitude')
# pl.ylabel(r'$\chi^2$')
# pl.show()

# pl.semilogy(tri['galb'][I], chi2, 'o')
# pl.xlabel('gal. latitude')
# pl.ylabel(r'$\chi^2$')
# pl.show()

# pl.loglog(d[I], chi2, 'o')
# pl.xlabel('d (pc)')
# pl.ylabel(r'$\chi^2$')
# pl.show()
