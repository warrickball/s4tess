#!/usr/bin/env python3

import numpy as np
from argparse import ArgumentParser

# parser = ArgumentParser()
# args = parser.parse_args()

tri = np.hstack([np.load('data/tri_v0_%02i.npy' %i) for i in range(13)])
atl = np.hstack([np.load('data/atl_v0_%02i.npy' %i) for i in range(13)])
ks = atl.dtype.names
X = atl.astype([(k, '>f8') for k in ks]).view('>f8').reshape((-1, len(ks)))
X[:,-2] = 1.0 - X[:,-2]  # Pmix -> 1-Pmix so ascending order gives largest Pmix first
X = np.flip(X, axis=1)  # reverse order of columns so unique() puts top-ranked stars first
X_unique, data_to_unique, unique_to_data = np.unique(X, axis=0, return_index=True, return_inverse=True)

np.save('data/tri_unique.npy', np.hstack([tri[i] for i in data_to_unique]))
np.save('data/atl_unique.npy', np.hstack([atl[i] for i in data_to_unique]))

sectors = -np.ones((len(X_unique), 13), dtype=int)

for i, j in enumerate(unique_to_data):
    sector, rank = np.divmod(i, 1000)
    sectors[j, sector] = rank

np.save('data/sectors.npy', sectors)
