#!/usr/bin/env python3

import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser(description=
"""Takes a TRILEGAL NumPy binary and writes MESA input files to
subfolders in the specified location."""
)

parser.add_argument('trilegal', type=str,
                    help="TRILEGAL NumPy binary with target stars")
parser.add_argument('basename', type=str,
                    help="""base output folder for Python format statement
                    (e.g. sector00/{:04d})""")
parser.add_argument('--inlist', type=str, default='inlist_run',
                    help="filename for MESA inlist (default inlist_run)")
parser.add_argument('--data', type=str, default='atl_data.txt',
                    help="filename for ATL data as text (default atl_data.txt)")
parser.add_argument('--verbose', '-v', action='store_const', const=True,
                    default=False, help="show progress")
args = parser.parse_args()

# approximate formula from Rodrigues et al. (2017)
Zsun = 0.01756
Ysun = 0.26618
Xsun = 1.0 - Zsun - Ysun
Yp = 0.2485
dY_dZ = 1.0068337129840546697 # 1.007

tri = np.load(args.trilegal)

for i, row in enumerate(tri):
    if args.verbose:
        print('\rCreating star {:d} of {:d}...'.format(i+1, len(tri)), end='')
        
    t = 10.**row['logAge']/1e9
    M = row['m_ini']
    Z = 10.**row['M_H']*Zsun  # as described in paper
    Y = Yp + dY_dZ*Z
    # Z_X = 10.**row['M_H']*Zsun/Xsun
    # X = (1.0 - Yp)/(1.0 + (1.0 + dY_dZ)*Z_X)
    # Z = Z_X*X
    # Y = 1.0 - X - Z
    with open(args.basename.format(i) + '/' + args.inlist, 'w') as f:
        f.writelines(['&star_job\n',
                      '/\n\n',
                      '&controls\n',
                      '    max_age = {:.16g}d9\n'.format(t),
                      '    initial_mass = {:.16g}d0\n'.format(M),
                      '    initial_Y = {:.16g}d0\n'.format(Y),
                      '    initial_Z = {:.16g}d0\n'.format(Z),
                      '    Zbase = {:.16g}d0\n'.format(Z),
                      '/\n\n',
                      '&pgstar\n',
                      '/\n'])

    with open(args.basename.format(i) + '/' + args.data, 'w') as f:
        for key in row.dtype.names:
            # print('{:s} = {:.16g}\n'.format(key, row[key]))
            f.write('{:>16s} = {:.16g}\n'.format(key, row[key]))

print()
