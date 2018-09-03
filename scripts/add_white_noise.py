#!/usr/bin/env python

import numpy as np
from tools import save_txt, load_txt
from argparse import ArgumentParser

parser = ArgumentParser(description=
"""Add white noise to a lightcurve using white noise level in metadata
file.""")

parser.add_argument('folder', type=str, help="folder containing files. "
                    "e.g. 'south/00000'")
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('--fmt', type=str, default='%16.7f',
                    help="format for output lightcurve (default='%16.7f')")
args = parser.parse_args()

def vprint(msg):
    if args.verbose:
        print(msg, end='', flush=True)

basename = args.folder + '/' + args.folder.split('/')[-1]

vprint('Loading metadata from %s/%s.meta... ' % (args.folder, basename))
meta = load_txt('%s.meta' % basename)
vprint('Done.\nLoading lightcurve from %s/%s.asc... ' % (args.folder, basename))

# asc = np.loadtxt('%s.asc' % basename)
with open('%s.asc' % basename, 'r') as f:
    asc = np.array([float(line) for line in f.readlines()])

vprint('Done.\nAdding white noise to lightcurve... ')
asc = asc + np.random.randn(len(asc))*meta['sigma']
vprint('Done.\nSaving lightcurve to %s/%s... ' % (args.folder, basename))

# np.savetxt('%s_WN.asc' % basename, asc, fmt='%16.7f')
with open('%s_WN.asc' % basename, 'w') as f:
    f.writelines([(args.fmt + '\n') % asci for asci in asc])

vprint('Done.\n')
