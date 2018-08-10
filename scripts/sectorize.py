#!/usr/bin/env python3

import numpy as np
from tools import sector_starts
from argparse import ArgumentParser

def vprint(*print_args, **kwargs):
    if args.verbose:
        print(*print_args, **kwargs)

parser = ArgumentParser()
parser.add_argument('filename', type=str,
                    help="filename containing truth table of which "
                    "stars occur in which sector")
parser.add_argument('--start', type=int, default=None,
                    help="Python-style index at which to start (default=first)")
parser.add_argument('--end', type=int, default=None,
                    help="Python-style index at which to end (default=last)")
# parser.add_argument('hemisphere', type=str,
#                     choices=['north', 'south'],
#                     help="which hemisphere ('north' or 'south')")
parser.add_argument('-v', '--verbose', action='store_const',
                    const=True, default=False)
args = parser.parse_args()

if 'north' in args.filename.lower():
    hemisphere = 'north'
    hemisector = 13
elif 'south' in args.filename.lower():
    hemisphere = 'south'
    hemisector = 0
else:
    raise ValueError("Couldn't work out which hemisphere applies. "
                     "Couldn't find 'north' or 'south' in filename.")

sectors = np.load(args.filename)
# N = 720*137//5  # (cadences/day)*(days/sector)

if args.start:
    if args.end:
        sectors = sectors[args.start:args.end]
    else:
        sectors = sectors[args.start:]
elif args.end:
    sectors = sectors[:args.end]

vprint('')

for star, row in enumerate(sectors):
    vprint('\rProcessing star %i...' % star, end='')

    try:
        with open('%s/%05i/%05i_WN.asc' % (hemisphere, star, star), 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        continue
        
    for sector, rank in enumerate(row):
        if rank < 0:
            continue
        else:
            start, end = sector_starts[hemisector+sector:hemisector+sector+2]
                
            with open('%s/%05i/%05i_WN_%02i_%04i.asc'
                      % (hemisphere, star, star, hemisector+sector, rank), 'w') as f:
                f.writelines(lines[start:end])
                
vprint('\nDone.')
