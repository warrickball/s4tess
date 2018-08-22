#!/usr/bin/env python3

import numpy as np
from tools import save_txt
from argparse import ArgumentParser

parser = ArgumentParser(description=
"""Loads the metadata in the specified file and writes it to text
files in the target folders."""
)
parser.add_argument('meta', type=str, help="filename of NumPy binary containing metadata")
parser.add_argument('filename', type=str,
                    help="format string for file, # will be replaced "
                    "by star number (e.g. #/#.meta)")
parser.add_argument('-N', type=int, default=[], nargs='+',
                    help="number of stars for which to create input "
                    "(default=all of them); if one number, do that many "
                    "stars starting from 0; if two numbers go from first "
                    "number to second *inclusive* (to match `seq`)")
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('--fmt', type=str, default='{:05d}', 
                    help="format statement for star number (default={:05d})")
parser.add_argument('--placeholder', type=str, default='#',
                    help="string to replace with star number in `filename`")
args = parser.parse_args()

def vprint(*print_args):
    if args.verbose:
        print(*print_args, end='', flush=True)

vprint('Loading metadata from %s... ' % args.meta)
meta = np.load(args.meta)
vprint('Done.\n')

try:
    end = args.N[1]+1
    start = args.N[0]
except IndexError:
    start = 0
    try:
        end = args.N[0]
    except IndexError:
        end = len(meta)

for i, row in enumerate(meta[start:end]):
    vprint('\rProcessing star %i of %i... ' % (i+1, end-start))
    save_txt(args.filename.replace(args.placeholder, args.fmt.format(i)), np.array(row))

vprint('Done.\n')
