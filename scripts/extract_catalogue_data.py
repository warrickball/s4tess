#!/usr/bin/env python3

import numpy as np
from argparse import ArgumentParser

def vprint(*print_args, **kwargs):
    "`print` only if `args.verbose` is True"
    if args.verbose:
        print(*print_args, **kwargs)

parser = ArgumentParser(description=
"""Takes a catalogue in a NumPy binary and writes data in plain text
to subfolders in the specified location."""
)

parser.add_argument('cat', type=str,
                    help="NumPy binary with catalogue of target stars")
parser.add_argument('output', type=str,
                    help="""base output filename for Python format statement
                    (e.g. sector00/{:04d}/cat_data.txt)""")
parser.add_argument('-N', type=int, default=[], nargs='+',
                    help="number of stars for which to create input "
                    "(default=all of them); if one number, do that many "
                    "stars starting from 0; if two numbers go from first "
                    "number to second *inclusive* (to match `seq`)")
parser.add_argument('--verbose', '-v', action='store_true',
                    help="show progress")
args = parser.parse_args()

cat = np.load(args.cat)

try:
    end = args.N[1]+1
    start = args.N[0]
except IndexError:
    end = len(cat)
    try:
        start = args.N[0]
    except IndexError:
        start = 0

vprint("Selecting rows from {:d} to {:d}...".format(start, end-1))

for i, row in enumerate(cat[start:end]):
    vprint('\rExtracting star {:d} of {:d}...'.format(i+1, end-start), end='')
    
    with open(args.output.format(start+i), 'w') as f:
        for key in row.dtype.names:
            f.write('{:>16s} = {:.16g}\n'.format(key, row[key]))

vprint('\nDone.')
