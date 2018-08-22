#!/usr/bin/env python3

import numpy as np
from tools import load_txt
from argparse import ArgumentParser
from scipy.interpolate import LinearNDInterpolator

def vprint(*print_args):
    if args.verbose:
        print(*print_args, end='', flush=True)

parser = ArgumentParser(description=
"""Creates MESA inlist files in the same folder as the specified
metadata file."""
)

parser.add_argument('meta', type=str,
                    help="text file containing metadata for target star")
parser.add_argument('--inlist', type=str, default='inlist_run',
                    help="filename for MESA inlist (default inlist_run)")
parser.add_argument('--verbose', '-v', action='store_true',
                    help="show progress")
parser.add_argument('--Tc', type=str,
                    default='/home/ADF/ballwh/work/s4tess/data/Tc.dat',
                    help="filename of central temperature data")
parser.add_argument('--nofit', action='store_true',
                    help="don't try to fit Teff and logL, just evolve to age")
args = parser.parse_args()

folder = '/'.join(args.meta.split('/')[:-1])

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

meta = load_txt(args.meta)

vprint("Processing... ")

Teff = 1e12 if args.nofit else 10.**meta['logTe']
log_L = 1e12 if args.nofit else meta['logL']
        
t = 10.**meta['logAge']/1e9
if not args.nofit:
    t = np.maximum(1.1*t, t + 0.2)

M = meta['m_ini']
Z = 10.**meta['M_H']*Zsun  # as described in paper
Y = Yp + dY_dZ*Z
# Z_X = 10.**meta['M_H']*Zsun/Xsun
# X = (1.0 - Yp)/(1.0 + (1.0 + dY_dZ)*Z_X)
# Z = Z_X*X
# Y = 1.0 - X - Z

Tc5 = Tc5_interpolator(M,Z)[()]

with open(folder + '/' + args.inlist, 'w') as f:
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

vprint('Done.\n')
