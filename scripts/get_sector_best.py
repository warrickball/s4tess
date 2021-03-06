#!/usr/bin/env python3

import numpy as np
from tools import tess_sectors
from argparse import ArgumentParser

def tess_fields_vector(lon, lat):
    Nsectors = 13

    n = 360.0/Nsectors
    try:
        sectors = np.full((len(lon), Nsectors), False)
    except TypeError:
        sectors = np.full((1, Nsectors), False)
        
    lon_range = np.minimum(12.0/abs(np.cos(np.radians(lat))), 180.0)

    for i in range(Nsectors):
        a = n*i
        d1 = abs(lon-a)
        d1 = np.minimum(d1, 360.0-d1)
        d2 = abs(lon - (a + 180.0)%360.0)
        d2 = np.minimum(d2, 360.0-d2)
        sectors[:,i] = ((d1 <= lon_range) & (6.0 <= lat)) | ((d2 <= lon_range) & (78.0 <= lat))
        # sectors[:,i] = (d1 <= lon_range and 6.0 <= lat) or (d2 <= lon_range and 78.0 <= lat)

    return sectors

parser = ArgumentParser(description=
"""Separates the stars in ATL file into sectors until all 26 sectors
have `N` targets in them.  Because the ATL file is sorted in
decreasing likelihood, these will be the `N` best targets."""
)

parser.add_argument('source', type=str, help="filename for input data")
parser.add_argument('output', type=str, help="base filename for output")
parser.add_argument('-N', '--Ntargets', type=int, default=1000,
                    help="number of targets per sector (default=1000)")
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('--start-lon', type=float, default=315.8,
                    # help="ecliptic longitude of first observing sector, "
                    # "in degrees (default=315.8)")
                    help="ecliptic longitude offset of the first observing sector, "
                    "in degrees (default=0)")
args = parser.parse_args()

def vprint(*print_args, **print_kwargs):
    if args.verbose:
        print(*print_args, **print_kwargs)

Nsectors = 26  # 13 North, 13 South

I = [[] for i in range(Nsectors)]

vprint('Loading data from %s...' % args.source)
data = np.load(args.source)

vprint('Selecting best stars in each sector...')
    
for i, row in enumerate(data):
    sectors_N = np.hstack([
        tess_fields_vector(row['ELon'] - args.start_lon, row['ELat'])[0],
        tess_fields_vector(row['ELon'] - args.start_lon, -row['ELat'])[0]
    ])
    
    for j, sector in enumerate(sectors_N):
        if sector and len(I[j]) < args.Ntargets:
            I[j].append(i)

    if np.all([len(i) >= args.Ntargets for i in I]): break

vprint('Saving output...')
        
for j, i in enumerate(I):
    np.save(args.output.format(j), data[i])

vprint('Done.')
