#!/usr/bin/env python3

import numpy as np
from numpy.lib.recfunctions import merge_arrays
from tools import tess_sectors, mean_v, sigma_v
from tools import sigma_schofield, sigma_stassun
from argparse import ArgumentParser

parser = ArgumentParser(description=
"""Separates the stars in an ATL file until all 26 sectors have at
least ``N`` targets.  Simultaneously cross-matches with TRILEGAL data
and computes additional metadata.  Because ATL is sorted by rank, each
sector will have its ``N`` best targets but many are duplicated in
multiple sectors, so the total number of targets is somewhat less than
26*N."""
)
parser.add_argument('atl', type=str, help="filename for ATL data")
parser.add_argument('tri', type=str, help="filename for TRILEGAL data")
parser.add_argument('meta', type=str, help="filename for output metadata table")
parser.add_argument('-N', type=int, default=1000,
                    help="number of targets per sector (default=1000)")
parser.add_argument('--noise-model', type=str, default='schofield',
                    choices=['schofield', 'stassun'],
                    help="choice of white noise model (default=schofield)")
parser.add_argument('-v', '--verbose', action='store_true')
args = parser.parse_args()

def vprint(*print_args):
    if args.verbose:
        print(*print_args, end='', flush=True)

vprint('Loading ATL data... ')
atl = np.load(args.atl)
atl = atl[np.abs(atl['ELat']) > 5.0]
vprint('Done.\n')

vprint('Loading TRILEGAL data... ')
tri = np.load(args.tri)
# try to make loop a bit faster by removing TRILEGAL targets that
# aren't within parameter space of ATL
tri = tri[(10.**tri['logL'] < np.max(atl['Lum'])*1.01) &
          (10.**tri['logL'] > np.min(atl['Lum'])*0.99) &
          (10.**tri['logR'] < np.max(atl['rad'])*1.01) &
          (10.**tri['logR'] > np.min(atl['rad'])*0.99) &
          (10.**tri['logTe'] < np.max(atl['teff']) + 10.0) &
          (10.**tri['logTe'] > np.min(atl['teff']) - 10.0)]
vprint('Done.\n')

ranks = -np.ones((26*args.N + 1000, 26), dtype=int)
J = -np.ones(len(ranks), dtype=int)

vprint('Computing sectors and cross-matching with TRILEGAL...\n')
vprint(''.join(['{:>12s}']*7).format('ATL l', 'TRI l',
                                     'ATL b', 'TRI b',
                                     'ATL Teff', 'TRI Teff', 'chi2'))
vprint('  min(N)\n')
for i, row in enumerate(atl):
    # vprint('\rProcessing star %i...' % (i+1))
    sectors = tess_sectors(row['ELon'], row['ELat'])
    
    if not np.any(sectors):
        continue

    # for j, sector in enumerate(sectors):
    #     if sector:
    #         ranks[i, j] = np.max(ranks[:,j]) + 1
    ranks[i, sectors] = np.max(ranks, axis=0)[sectors] + 1

    allchi2 = (row['GLon']-tri['gall'])**2 + (row['GLat']-tri['galb'])**2

    J[i] = np.argmin(allchi2)
    vprint('{:>12.6f}{:>12.6f}{:>12.6f}{:>12.6f}{:>12.3f}{:>12.3f}'.format(
        row['GLon'], tri['gall'][J[i]], row['GLat'], tri['galb'][J[i]],
        row['teff'], 10.**tri[J[i]]['logTe']))

    vprint('%12.3e' % np.min(allchi2))

    minmaxrank = np.min(np.max(ranks, axis=0))
    vprint('%8i\n' % minmaxrank)
    
    if np.min(allchi2) > 1e-8:
        print('ERROR: allchi2 > 1e-8')
        exit
    
    if minmaxrank >= args.N:
        break

vprint('Done.\n')

I = np.where(np.any(ranks >= 0, axis=1))[0]

vprint('White noise model is `%s`.' % args.noise_model)
vprint('Combining datasets... ')
dtype = [('gall', float), ('galb', float), ('Gc', int), ('logAge', float),
         ('logTe', float), ('logL', float),
         ('M_H', float), ('m_ini', float), ('mu0', float), ('Av', float),
         ('comp', int), ('Mass', float), ('imag', float), ('TESSmag', float),
         ('ELon', float), ('ELat', float),
         ('Pdet_fixedBeta', float), ('SNR_fixedBeta', float),
         ('Pdet_varyBeta', float), ('SNR_varyBeta', float), ('P_mix', float),
         ('Rank_Pmix', float), ('vr', float), ('sigma', float)] + [('rank_%02i' % i, int) for i in range(26)]

meta = np.zeros(len(I), dtype=dtype)

for k in meta.dtype.names:
    if k in atl.dtype.names:
        meta[k] = atl[I][k]
    elif k in tri.dtype.names:
        meta[k] = tri[J[I]][k]
    elif k == 'vr':
        lon = tri[J[I]]['gall']
        meta[k] = mean_v(lon) + np.random.randn()*sigma_v(lon)
    elif k == 'sigma':
        if args.noise_model == 'stassun':
            meta[k] = sigma_stassun(atl[I]['imag'])
        elif args.noise_model == 'schofield':
            meta[k] = sigma_schofield(tri[J[I]]['imag'], atl[I]['teff'], 
                                      atl[I]['ELon'], atl[I]['ELat'], 
                                      atl[I]['GLon'], atl[I]['GLat'])
    elif k.startswith('rank_'):
        continue
    else:
        raise ValueError('no field of name %s' % k)

for i in range(26):
    meta['rank_%02i' % i] = ranks[I][:,i]

vprint('Done.\nSaving output to %s... ' % args.meta)
np.save(args.meta, meta)
vprint('Done.\n')
