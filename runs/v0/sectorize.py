#!/usr/bin/env python3

import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_const',
                    const=True, default=False)
args = parser.parse_args()

sectors = np.load('data/sectors.npy')
N = 19440

if args.verbose:
    print('')

for star, row in enumerate(sectors):
    if args.verbose:
        print('\rProcessing star %i...' % star, end='')
    try:
        with open('unique/%05i/unique_%05i_WN.asc' % (star, star), 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        continue
        
    for sector, rank in enumerate(row):
        if rank < 0:
            continue
        else:
            with open('unique/%05i/unique_%05i_WN_%02i_%04i.asc'
                      % (star, star, sector, rank), 'w') as f:
                f.writelines(lines[sector*N:sector*N+N])
                
if args.verbose:
    print('\nDone.')
