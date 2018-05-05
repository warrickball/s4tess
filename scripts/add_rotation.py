#!/usr/bin/env python

import numpy as np
from scipy.optimize import fsolve
from tomso import gyre, mesa
from argparse import ArgumentParser

def BV_to_Teff(BV):
    # coefficients as in Table 2, which is revers of polyval input
    p = [3.979145106714099, -0.654992268598245, 1.740690042385095,
         -4.608815154057166, 6.792599779944473, -5.396909891322525,
         2.192970376522490, -0.359495739295671]
    return 10.**np.polyval(p[::-1], BV)

def Teff_to_BV(Teff):
    return fsolve(lambda z: BV_to_Teff(z)-Teff, np.zeros_like(Teff))

def P_mamabrand(Teff, t, noisy=False):
    # Mamajek & Hillebrand (2008), eq. (12)--(14), period in days
    # http://adsabs.harvard.edu/abs/2008ApJ...687.1264M
    k = 1.0 if noisy else 0.0
    a = 0.407 + 0.021*noisy*np.random.randn()
    b = 0.325 + 0.024*noisy*np.random.randn()
    c = 0.495 + 0.010*noisy*np.random.randn()
    n = 0.566 + 0.008*noisy*np.random.randn()
    return a*(Teff_to_BV(Teff)-c)**b*t**n

parser = ArgumentParser(description="""
Add rotation to the specified model.""")
parser.add_argument('model', type=str,
                  help="GYRE model to which to add rotation.")
parser.add_argument('history', type=str,
                  help="MESA history (for age and other properties)")
parser.add_argument('-o', '--output', type=str, default=None,
                 help="Filename for output GYRE model "
                 "(default=append .rot to source model filename)")
parser.add_argument('--no-core-noise', action='store_const',
                    const=True, default=False,
                    help="don't randomly sample uncertainty in core rotation rate")
parser.add_argument('--no-env-noise', action='store_const',
                    const=True, default=False,
                    help="don't randomly sample uncertainty in envelope rotation rate")
args = parser.parse_args()

gyre_header, profile = gyre.load_gyre(args.model)
mesa_header, history = mesa.load_history(args.history)

P = P_mamabrand(10.**history['log_Teff'][-1], history['star_age'][-1]/1e6,
                noisy=not args.no_env_noise)*86400.0
Omega_env = 2.*np.pi/P
profile['omega'] = Omega_env

if history['center_h1'][-1] < 1e-4 and np.log10(history['gravity'][-1]) < 3.8:
    # fit to Mosser et al. (2012) data for RGB stars, Dnu > 12 uHz
    if args.no_core_noise:
        Omega_core = 0.375
    else:
        Omega_core = 0.375 + 0.1045*np.random.randn()
        
    Omega_core = Omega_core*1e-6*2.*np.pi  # nHz -> radians/s
    Omega_core = max(Omega_core, Omega_env)  # force Omega_core > Omega_env
    
    # find the core/envelope boundary
    k_core = [row['k'] for row in profile
              if row['r']/gyre_header['R'] < 0.98 and row['N2'] > 0][-1]
    profile['omega'][:k_core] = Omega_core

if args.output is None:
    output = args.model + '.rot'
else:
    output = args.output

gyre.save_gyre(output, gyre_header, profile)
