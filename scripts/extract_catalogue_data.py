#!/usr/bin/env python3

import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser(description=
"""Takes a catalogue in a NumPy binary and writes data in plain text
to subfolders in the specified location."""
)

parser.add_argument('cat', type=str,
                    help="NumPy binary with catalogue of target stars")
parser.add_argument('output', type=str,
                    help="""base output filename for Python format statement
                    (e.g. sector00/{:04d}/cat_data.txt)""")
parser.add_argument('--verbose', '-v', action='store_const', const=True,
                    default=False, help="show progress")
args = parser.parse_args()

cat = np.load(args.cat)

for i, row in enumerate(cat):
    if args.verbose:
        print('\rCreating star {:d} of {:d}...'.format(i+1, len(cat)), end='')
    
    with open(args.output.format(i), 'w') as f:
        for key in row.dtype.names:
            f.write('{:>16s} = {:.16g}\n'.format(key, row[key]))

if args.verbose:
    print('\nFinished')
