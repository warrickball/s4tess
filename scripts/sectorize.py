#!/usr/bin/env python3

import numpy as np
from tools import load_txt, sector_starts
from argparse import ArgumentParser

def vprint(*print_args, **kwargs):
    if args.verbose:
        print(*print_args, **kwargs)

parser = ArgumentParser()
parser.add_argument('folder', type=str,
                    help="filename containing truth table of which "
                    "stars occur in which sector")
parser.add_argument('-v', '--verbose', action='store_const',
                    const=True, default=False)
args = parser.parse_args()

basename = args.folder.split('/')[-1]
meta = load_txt(args.folder + '/' + basename + '.meta')
ranks = [meta['rank_%02i' % i] for i in range(26)]

vprint('Separating lightcurve in %s in sectors... ' % args.folder)

with open(args.folder + '/' + basename + '_WN.asc', 'r') as f:
    lines = f.readlines()

for sector, rank in enumerate(ranks):
    if rank < 0:
        continue
    else:
        start, end = sector_starts[sector:sector+2]
        with open('%s/%s_WN_%02i_%04i.asc'
                  % (args.folder, basename, sector, rank), 'w') as f:
            f.writelines(lines[start:end])
                
vprint('\nDone.')
