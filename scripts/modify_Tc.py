#!/usr/bin/env python3

from argparse import ArgumentParser

parser = ArgumentParser(description=
"""Small script to adjust ``pre_ms_T_c`` if initial model failed to
converge.""")
parser.add_argument('inlist', type=str, help="filename of inlist")
parser.add_argument('T_c', type=float, 
                    help="initial central temperature (pre_ms_T_c in MESA)")
parser.add_argument('-d', action='store_const', const=True, default=False, 
                    help="input value is a temperature difference rather than a temperature")
parser.add_argument('-o', '--output', type=str, default=None,
                    help="filename for output (defaults to input inlist)")
args = parser.parse_args()

with open(args.inlist, 'r') as f:
    lines = [line.split('=') for line in f.readlines()]

i = [i for i, line in enumerate(lines) if 'pre_ms_T_c' in line[0]][0]

if args.d:
    T_c = float(lines[i][1].replace('d','e')) + args.T_c
else:
    T_c = args.T_c

if T_c >= 1e6:
    raise ValueError('central temperature {:.3e} K, must be less than 1e6 K'.format(T_c))
    
lines[i][1] = ' {:.16g}d5\n'.format(T_c/1e5)

output = args.inlist if args.output is None else args.output
with open(output, 'w') as f:
    f.writelines(['='.join(line) for line in lines])
