#!/usr/bin/env python

import numpy as np
from tomso import gyre, io
from argparse import ArgumentParser
from solioak import scaling
import AADG3

parser = ArgumentParser()
parser.add_argument('summary', type=str, help='input GYRE summary')
parser.add_argument('namein', type=str, help='name of output Fortran namelist')
parser.add_argument('namecon', type=str, help='name of output .con file')
parser.add_argument('namerot', type=str, help='name of output .rot file')
parser.add_argument('--splitting', type=float, default=0.0,
                    help='constant rotation splitting')
args = parser.parse_args()

header, summary = gyre.load_summary(args.summary)

Lsun = 3.844e33
Rsun = 6.96568e10
Msun = 1.989e33
Teff_sun = 5777.0
numax_sun = 3090.0
Dnu_sun = 135.1
sig_sun = 60.0

M = header['M_star']/Msun  # 1.0  # Msun
R = header['R_star']/Rsun  # 1.0  # Rsun
L = header['L_star']/Lsun  # 1.0  # Lsun
Teff = (L/R**2)**0.25*Teff_sun
numax = numax_sun*(M/R**2/(Teff/Teff_sun)**0.5)
Dnu = Dnu_sun*np.sqrt(M/R**3)
sig = np.sqrt(L**2/M**3/(Teff/Teff_sun)**5.5*(numax/numax_sun))*sig_sun

nml = {'user_seed': 0, 'cadence': 120.0, 'n_cadences': 19440, 'n_relax': 4320, 'n_fine': 50,
       'sig': 60.0, 'rho': 0.45, 'tau': 250.0, 'inclination': 90.0,
       'pcyc': 100.0, 'phi': 0.0, 'nuac': 0.0,
       'p(1)': 1.52355, 'p(2)': 0.565349, 'p(3)': 0.0361707,
       'ass_init': True, 'namecon': args.namecon, 'namerot': args.namerot,
       'nameout': args.namecon.replace('.con', '.asc')}

nml['tau'] = 250.0/(numax/3090.)
nml['sig'] = sig

AADG3.save_nml(args.namein, nml)

l = summary['l'].astype(int)
n = summary['n_pg'].astype(int)
nu = summary['Refreq']
E = summary['E_norm']

# all As are A_rms
# 
Amax_sun = 2.1        # ppm
Tred_sun = 8907.0
dT = 1250.0
    
Tred = Tred_sun*L**-0.093
beta = 1.0 - np.exp((Teff-Tred)/dT)
Amax = Amax_sun*beta*L/M*(Teff/Teff_sun)**-2
Henv = Amax**2/Dnu
wenv = 0.66*numax**0.88
if Teff > Teff_sun:
    wenv *= 1.0 + 6e-4*(Teff-Teff_sun)
    
cenv = wenv/2./np.sqrt(2.*np.log(2.))

Q = E/np.interp(nu, nu[l==0], E[l==0])

width = scaling.lund_width(nu, numax)/Q

H = Henv*np.exp(-(nu-numax)**2/2./cenv**2)/Q
amp2 = H*Dnu

np.savetxt(args.namecon, np.vstack([l, n, nu, width, amp2, 0.0*nu]).T,
           fmt=['%8i','  %7i','  %12.3e','  %12.7e','  %12.7e','  %12.8e'])

with open(args.namerot, 'w') as f:
    for ni, li in zip(n, l):
        for m in range(1, li+1):
            f.write('%8i%8i%8i%9.3f\n' % (ni, li, m, args.splitting))
    
