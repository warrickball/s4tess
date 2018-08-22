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
args = parser.parse_args()

def vprint(msg):
    if args.verbose:
        print(msg, end='', flush=True)

basename = args.folder + '/' + args.folder.split('/')[-1]

vprint('Loading metadata... ')
meta = load_txt('%s.meta' % basename)
vprint('Done.\nAdding white noise to lightcurve... ')
asc = np.loadtxt('%s.asc' % basename)
asc = asc + np.random.randn(len(asc))*meta['sigma']
np.savetxt('%s_WN.asc' % basename, asc, fmt='%16.7f')
vprint('Done.\n')
