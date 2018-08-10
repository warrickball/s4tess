#!/usr/bin/env python3

import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser(description=
"""Takes lists of top targets in each sector and saves a NumPy array
of unique targets along with a NumPy array specifying in which sectors
those targets occur.""")
parser.add_argument('hemisphere', type=str,
                    choices=['north', 'south'],
                    help="which hemisphere ('north' or 'south')")
args = parser.parse_args()

if args.hemisphere is 'north':
    sector_range = range(13)
elif args.hemisphere is 'south':
    sector_range = range(13, 26)
else:
    raise ValueError('invalid choice for args.hemisphere\n(this should never be evaluated...)')

tri = np.hstack([np.load('data/tri_v0_%02i.npy' %i) for i in sector_range])
atl = np.hstack([np.load('data/atl_v0_%02i.npy' %i) for i in sector_range])
ks = atl.dtype.names
X = atl.astype([(k, '>f8') for k in ks]).view('>f8').reshape((-1, len(ks)))
X[:,-2] = 1.0 - X[:,-2]  # Pmix -> 1-Pmix so ascending order gives largest Pmix first
X = np.flip(X, axis=1)  # reverse order of columns so unique() puts top-ranked stars first
X_unique, data_to_unique, unique_to_data = np.unique(X, axis=0, return_index=True, return_inverse=True)

np.save('data/tri_%s.npy' % args.hemisphere, np.hstack([tri[i] for i in data_to_unique]))
np.save('data/atl_%s.npy' % args.hemisphere, np.hstack([atl[i] for i in data_to_unique]))

sectors = -np.ones((len(X_unique), 13), dtype=int)

for i, j in enumerate(unique_to_data):
    sector, rank = np.divmod(i, 1000)
    sectors[j, sector] = rank

np.save('data/sectors_%s.npy' % args.hemisphere, sectors)
