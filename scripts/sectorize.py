#!/usr/bin/env python3

import numpy as np
from tools import load_txt, sector_starts
from argparse import ArgumentParser

def vprint(*print_args):
    if args.verbose:
        print(*print_args, end='', flush=True)

parser = ArgumentParser()
parser.add_argument('folder', type=str,
                    help="folder containing output from AADG3")
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('--WN', action='store_true',
                    help="use the lightcurve with white noise added "
                    "(i.e. with `WN` in filename)")
args = parser.parse_args()

basename = args.folder.split('/')[-1]
meta = load_txt(args.folder + '/' + basename + '.meta')
ranks = [meta['rank_%02i' % i] for i in range(26)]

vprint('Separating lightcurve in %s in sectors... ' % args.folder)

suffix = '_WN' if args.WN else ''

with open(args.folder + '/' + basename + suffix + '.asc', 'r') as f:
    lines = f.readlines()

for sector, rank in enumerate(ranks):
    if rank < 0:
        continue
    else:
        hemi = sector//13*13
        start, end = sector_starts[sector:sector+2]
        start -= sector_starts[hemi]
        end -= sector_starts[hemi]
        with open('%s/%s%s_%02i_%04i.asc'
                  % (args.folder, basename, suffix, sector, rank), 'w') as f:
            f.writelines(lines[start:end])
                
vprint('Done.\n')
