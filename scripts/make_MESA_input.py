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
                    help="""base output filename for Python format statement
                    (e.g. sector00/{:4d}/inlist_run)""")
args = parser.parse_args()

# approximate formula from Rodrigues et al. (2017)
Zsun = 0.01756
Yp = 0.2485
dY_dZ = 1.007

tri = np.load(args.trilegal)

for i, row in enumerate(tri):
    t = 10.**row['logAge']/1e9
    M = row['m_ini']
    Z = 10.**row['M_H']*Zsun
    Y = Yp + dY_dZ*Z
    with open(args.basename.format(i), 'w') as f:
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
