#!/usr/bin/env python3

import numpy as np
from argparse import ArgumentParser
from scipy.interpolate import LinearNDInterpolator

def vprint(*print_args, **kwargs):
    "`print` only if `args.verbose` is True"
    if args.verbose:
        print(*print_args, **kwargs)

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
parser.add_argument('--tri-data', type=str, default=None,
                    help="save TRILEGAL data as text to this filename")
parser.add_argument('--verbose', '-v', action='store_true',
                    help="show progress")
parser.add_argument('--Tc', type=str,
                    default='/home/ADF/ballwh/work/s4tess/data/Tc.dat',
                    help="filename of central temperature data")
parser.add_argument('--nofit', action='store_true',
                    help="don't try to fit Teff and logL, just evolve to age")
parser.add_argument('-N', type=int, default=[], nargs='+',
                    help="number of stars for which to create input "
                    "(default=all of them); if one number, do that many "
                    "stars starting from 0; if two numbers go from first "
                    "number to second *inclusive* (to match `seq`)")
args = parser.parse_args()

# get central temperature data
Tc = np.genfromtxt(args.Tc, names=True)
Tc5_interpolator = LinearNDInterpolator(np.vstack([Tc['M'], Tc['Z']]).T,
                                        np.around(Tc['Tc']/1e5))

# approximate formula from Rodrigues et al. (2017)
Zsun = 0.01756
Ysun = 0.26618
Xsun = 1.0 - Zsun - Ysun
Yp = 0.2485
dY_dZ = 1.0068337129840546697 # 1.007

tri = np.load(args.trilegal)

try:
    end = args.N[1]+1
    start = args.N[0]
except IndexError:
    start = 0
    try:
        end = args.N[0]
    except IndexError:
        end = len(tri)

vprint("Selecting rows from {:d} to {:d}...".format(start, end-1))
        
for i, row in enumerate(tri[start:end]):
    vprint('\rCreating star {:d} of {:d}...'.format(i+1, end-start), end='')

    Teff = 1e12 if args.nofit else 10.**row['logTe']
    log_L = 1e12 if args.nofit else row['logL']
        
    t = 10.**row['logAge']/1e9
    if not args.nofit:
        t = np.maximum(1.1*t, t + 0.2)

    M = row['m_ini']
    Z = 10.**row['M_H']*Zsun  # as described in paper
    Y = Yp + dY_dZ*Z
    # Z_X = 10.**row['M_H']*Zsun/Xsun
    # X = (1.0 - Yp)/(1.0 + (1.0 + dY_dZ)*Z_X)
    # Z = Z_X*X
    # Y = 1.0 - X - Z

    Tc5 = Tc5_interpolator(M,Z)[()]

    with open(args.basename.format(start+i) + '/' + args.inlist, 'w') as f:
        f.writelines(['&star_job\n',
                      '    pre_ms_T_c = {:.16g}d5\n'.format(Tc5),
                      '/\n\n',
                      '&controls\n',
                      '    max_age = {:.16g}d9\n'.format(t),
                      '    initial_mass = {:.16g}d0\n'.format(M),
                      '    initial_Y = {:.16g}d0\n'.format(Y),
                      '    initial_Z = {:.16g}d0\n'.format(Z),
                      '    Zbase = {:.16g}d0\n'.format(Z),
                      '\n',
                      '    x_ctrl(1) = {:.16g}d0\n'.format(Teff),
                      '    x_ctrl(2) = {:.16g}d0\n'.format(log_L),
                      '/\n\n',
                      '&pgstar\n',
                      '/\n'])

    if args.tri_data:
        with open(args.basename.format(start+i) + '/' + args.tri_data, 'w') as f:
            for key in row.dtype.names:
                # print('{:s} = {:.16g}\n'.format(key, row[key]))
                f.write('{:>16s} = {:.16g}\n'.format(key, row[key]))

print()
