#!/usr/bin/env python

import numpy as np
from tomso import gyre
from tools import save_txt, load_txt, sector_lengths
from argparse import ArgumentParser
# from solioak import scaling
from collections import OrderedDict
from hashlib import sha1  # overkill?
import AADG3

def logW_meta(nu, numax, Teff, a_alpha, b_alpha, c_alpha,
              a_Walpha, b_Walpha, c_Walpha,
              a_DWdip, b_DWdip, c_DWdip,
              a_nudip, b_nudip, c_nudip,
              a_Wdip, b_Wdip, c_Wdip):
    alpha = a_alpha + b_alpha*Teff + c_alpha*numax
    Walpha = a_Walpha + b_Walpha*Teff + c_Walpha*numax
    DWdip = a_DWdip + b_DWdip*Teff + c_DWdip*numax
    nudip = a_nudip + b_nudip*Teff + c_nudip*numax
    Wdip = a_Wdip + b_Wdip*Teff + c_Wdip*numax

    output = alpha*np.log(nu/numax) + np.log(Walpha)
    output[DWdip < 1.0] += (np.log(DWdip)/(1. + (2*np.log(nu/nudip)/np.log(Wdip))**2))[DWdip < 1.0]
    return output


parser = ArgumentParser()
# parser.add_argument('summary', type=str, help='input GYRE summary')
# parser.add_argument('atl', type=str, help='input ATL data (for galactic longitude)')
# parser.add_argument('baseout', type=str, help='basename for output files')
parser.add_argument('folder', type=str, help="folder containing files to manipulate "
                    "(must exclude trailing /)")
parser.add_argument('--splitting', type=float, default=-1,
                    help="constant rotation splitting, ignored if < 0, "
                    "in which case rotation is taken from summary "
                    "if available, otherwise rotation = 0 (default -1)")
parser.add_argument('--min-H', type=float, default=-1,
                    help="minimum for height H (default=-1, i.e. keep all modes)")

parser.add_argument('-v', '--verbose', action='store_true')
args = parser.parse_args()

basename = args.folder.split('/')[-1]
meta = load_txt(args.folder + '/%s.meta' % basename)

header, summary = gyre.load_summary(args.folder + '/gyre_summary.txt')

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

modes_filename = basename + '.con'
rotation_filename = basename + '.rot'

np.random.seed(int(sha1(('azaza' + args.folder + 'bybyb').encode('utf-8')).hexdigest(), 16)%2**32)
warmup = np.random.randint(2**16, size=100)

namelist = OrderedDict()
namelist['user_seed'] = np.random.randint(100, 2**28-1)
namelist['cadence'] = 120.0
# namelist['n_cadences'] = 13*720*137//5  # n_sectors*(cadences/day)*(days/sector=27.4d)
# if 'south' in args.folder.lower():
if np.any([meta['rank_%02i' % i] > -1 for i in range(0, 13)]):
    namelist['n_cadences'] = sum(sector_lengths[0:13])
# elif 'north' in args.folder.lower():
elif np.any([meta['rank_%02i' % i] > -1 for i in range(13, 26)]):
    namelist['n_cadences'] = sum(sector_lengths[13:26])
else:
    raise ValueError("Can't tell if star is in northern or southern hemisphere! "
                     "Expected to find either 'north' or 'south' in folder name.")

namelist['n_relax'] = 6*720
namelist['n_fine'] = 50
namelist['sig'] = sig
namelist['rho'] = 0.45
namelist['tau'] = 250.0/(numax/3090.)
namelist['inclination'] = np.degrees(np.arccos(np.random.rand()))
namelist['cycle_period'] = 100.0
namelist['cycle_phase'] = 0.0
namelist['nuac'] = 0.0
namelist['p(1)'] = 1.52355
namelist['p(2)'] = 0.565349
namelist['p(3)'] = 0.0361707
namelist['add_granulation'] = True
namelist['modes_filename'] = modes_filename
namelist['rotation_filename'] = rotation_filename
namelist['output_filename'] = basename + '.asc'

AADG3.save_namelist(args.folder + '/' + basename + '.in', namelist)

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

# width = scaling.lund_width(nu, numax)/Q
# using 25 of Guy's red giants
p25 = [-3.71033159e+00,  1.07268220e-03,  1.88285544e-04, -7.20902433e+01,
        1.54336225e-02,  9.10050555e-04, -2.26620472e-01,  5.08279583e-05,
        2.71537654e-06, -2.18970560e+03,  4.30163078e-01,  8.42663954e-01,
       -5.63927133e-01,  1.13801669e-04,  1.31215914e-04]
width = np.exp(logW_meta(nu, numax, Teff, *p25))/Q

H = Henv*np.exp(-(nu-numax)**2/2./cenv**2)/Q
amp2 = H*Dnu
I = H > args.min_H

# Doppler shift the frequencies
TAU = 2*np.pi
def func(x, a0, a1, p1, a2, p2):
    return a0 + a1*np.sin(TAU*(x/360.0 + p1)) + a2*np.sin(TAU*(x/180.0 + p2))

v = meta['vr']
c = 299792.458

nu = nu*np.sqrt((1-v/c)/(1+v/c))

np.savetxt(args.folder + '/' + modes_filename, 
           np.vstack([l, n, nu, width, amp2, 0.0*nu]).T[I],
           fmt=['%2i','  %5i','  %12.7e','  %12.7e','  %12.7e','  %12.8e'])

if args.splitting >= 0:
    dnu_rot = args.splitting*np.ones_like(nu)
else:
    try:
        dnu_rot = summary['dfreq_rot']
    except ValueError:
        dnu_rot = np.zeros_like(nu)
        

with open(args.folder + '/' + rotation_filename, 'w') as f:
    for ni, li, dnu_roti in zip(n[I], l[I], dnu_rot[I]):
        for m in range(1, li+1):
            f.write('%5i%3i%3i%12.7f\n' % (ni, li, m, dnu_roti))

