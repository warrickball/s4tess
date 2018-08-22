#!/usr/bin/env python

import numpy as np
from tools import save_txt, load_txt
from argparse import ArgumentParser

parser = ArgumentParser(description=
"""Add white noise to a lightcurve, either using the formula in
Stassun et al. (2017) or Mat Schofield and Bill Chaplin's formula.

All files are named based on the following convention with the
`basename` argument:

* the ATL data is in `basename.atl`,
* the TRILEGAL data is in `basename.tri`,
* the white noise value is saved to `basename.xtras`,
* the input light curve (without white noise) is in `basename.asc`, and
* the output light curve (after adding white noise) is in `basename_WN.asc`.""")

parser.add_argument('folder', type=str, help="folder containing files. "
                    "e.g. 'south/00000'")
# parser.add_argument('atl', type=str,
#                     help="plain text file with ATL data")
# parser.add_argument('tri', type=str,
#                     help="plain text file with TRILEGAL data")
# parser.add_argument('lc_in', type=str,
#                     help='filename of input lightcurve')
# parser.add_argument('lc_out', type=str,
#                     help='filename of output lightcurve')
# parser.add_argument('--model', type=str, default='schofield',
#                     choices=['schofield', 'stassun'],
#                     help="choice of noise model")
parser.add_argument('-v', '--verbose', action='store_true')
args = parser.parse_args()

# def stassun(Imag):
#     # see Stassun et al. (2017)
#     # https://arxiv.org/abs/1706.00495v3

#     # if Imag is less than local minimum aroung 5,
#     # use flat line that roughly corresponds to systematic noise
#     Imag = np.maximum(Imag, 4.97756051)
#     F = 4.73508403525e-5
#     E = -0.0022308015894
#     D = 0.0395908321369
#     C = -0.285041632435
#     B = 0.850021465753
#     lnA = 3.29685004771
#     lnsigma = lnA + B*Imag + C*Imag**2 + D*Imag**3 + E*Imag**4 + F*Imag**5
#     return np.exp(lnsigma)


# def schofield(Imag, Teff, ELon=0.0, ELat=30.0, GLon=96.0, GLat=-30.0,
#               cadence=120.0, subcadence=2.0, npix_aper=4, frac_aper=0.76,
#               e_pix_ro=10, geom_area=60.0, pix_scale=21.1, sys_limit=0.0):
#     # from Mat Schofield's ATL code
#     # "from the TESS Wiki"...

#     omega_pix = pix_scale**2
#     n_exposures = cadence/subcadence

#     # electrons from the star
#     megaph_s_cm2_0mag = 1.6301336 + 0.14733937*(Teff-5000.0)/5000.0
#     e_star = 10.0**(-0.4*Imag)*1e6*megaph_s_cm2_0mag*geom_area*cadence*frac_aper
#     e_star_sub = e_star*subcadence/cadence

#     # e/pix from zodi
#     dlat = (abs(ELat) - 90.0)/90.0
#     vmag_zodi = 23.345 - 1.148*dlat**2
#     e_pix_zodi = 10.0**(-0.4*(vmag_zodi-22.8))*2.39e-3*geom_area*omega_pix*cadence

#     # e/pix from background stars
#     dlat = abs(GLat)/40.0

#     dlon = GLon
#     # q = np.where(dlon > 180.0)
#     # if len(q[0]) > 0:
#     # 	dlon[q] = 360.0 - dlon[q]
#     if dlon > 180.0:
#         dlon = 360.0 - dlon

#     dlon = abs(dlon)/180.0
#     p = [18.97338, 8.833, 4.007, 0.805]
#     Imag_bgstars = p[0] + p[1]*dlat + p[2]*dlon**p[3]
#     e_pix_bgstars = 10.0**(-0.4*Imag_bgstars)*1.7e6*geom_area*omega_pix*cadence

#     # compute noise sources
#     noise_star = np.sqrt(e_star)/e_star
#     noise_sky  = np.sqrt(npix_aper*(e_pix_zodi + e_pix_bgstars))/e_star
#     noise_ro   = np.sqrt(npix_aper*n_exposures)*e_pix_ro/e_star
#     noise_sys  = 0.0*noise_star + sys_limit/1e6/np.sqrt(cadence/3600.0)

#     noise = np.sqrt(noise_star**2 + noise_sky**2 + noise_ro**2 + noise_sys**2)
#     return noise*1e6*np.sqrt(cadence/3600.0)  # parts x sqrt(hr) -> ppm * sqrt(cadence)


def vprint(msg):
    if args.verbose:
        print(msg, end='', flush=True)

basename = args.folder + '/' + args.folder.split('/')[-1]

vprint('Loading metadata... ')
meta = load_txt('%s.meta' % basename)
vprint('Done.\n')

vprint('Done.\nAdding white noise to lightcurve... ')
    
asc = np.loadtxt('%s.asc' % basename)
asc = asc + np.random.randn(len(asc))*meta['sigma']
np.savetxt('%s_WN.asc' % basename, asc, fmt='%16.7f')
vprint('Done.\n')
