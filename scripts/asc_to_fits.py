#!/usr/bin/env python

import numpy as np
import AADG3
from tools import load_txt, sector_starts
from tomso import gyre, mesa
from astropy.io import fits
from datetime import datetime
from argparse import ArgumentParser

parser = ArgumentParser(description="""Takes results of S4TESS
simulation and writes a vaguely MAST-like FITS file.""")
parser.add_argument('asc', type=str,
                    help="Filename of output timeseries (usually *_WN_*.asc).")
parser.add_argument('-o', '--output', type=str, default=None,
                    help="Alternative filename for output FITS file.  "
                    "By default, the output filename is same as the input filename "
                    "with `.asc` replaced with `.fits`.  If the input filename "
                    "doesn't end with `.asc`, the output file is just the input "
                    "filename with `.fits` appended.")
parser.add_argument('--skip-test', action='store_true',
                    help="don't test the FITS file after it's written")
# args = parser.parse_args(['../runs/v0/south/00000/00000_WN_11_0000.asc'])
args = parser.parse_args()

forms = {'J': 'column format: signed 32-bit integer',
         'D': 'column format: 64-bit floating point',
         'E': 'column format: 32-bit floating point'}
cads_per_sector = 720*137//5  # (cadences/day)*(days/sector)

ascname = args.asc.split('/')[-1]
folder = '/'.join(args.asc.split('/')[:-1])
basename = folder + '/' + args.asc.split('/')[-2]
ID, sector, sec_rank = map(int, [ascname.split('.')[0].split('_')[i] for i in [0,2,3]])

if args.output:
    fitsname = args.output
else:
    if args.asc.endswith('.asc'):
        fitsname = args.asc[:-4] + '.fits'
    else:
        fitsname = args.asc + '.fits'

# tri = {}  # TRILEGAL data
# with open(basename + '.tri', 'r') as f:
#     for line in f.readlines():
#         k, v = line.split('=')
#         tri[k.strip()] = float(v)
meta = load_txt(basename + '.meta')

meta_comments = {
    'gall': 'galactic longitude',
    'galb': 'galactic latitude',
    'Gc': 'galactic component: 1--4 = thin disc, thick disc, halo, bulge',
    'logAge': '[yr] log10 age',
    'M_H': 'metallicity',
    'm_ini': '[solar masses] initial stellar mass',
    'mu0': 'distance modulus',
    'Av': 'reddening',
    'comp': 'primary or secondary (1 or 2) in a binary',
    'Mass': '[solar masses] current stellar mass',
    'logL': '[solar luminosities] log10 luminosity',
    'logTe': '[K] log10 effective temperature',
    'logg': '[cm/s2] surface gravity',
    'label': '',
    'logR': '[solar radii] log10 stellar radius',
    'logrho': '[g/cm3] log10 mean stellar density',
    'numax': '[uHz] frequency of maximum oscillation power',
    'ddeltanu': '',
    'deltanu': '[uHz] large separation',
    'pspaci': '',
    'pcoupl': '',
    'mbolmag': '',
    'TESSmag': 'TESS magnitude',
    'Jmag': '',
    'Hmag': '',
    'Ksmag': '',
    'Keplermag': 'Kepler magnitude',
    'gmag': 'g-band magnitude',
    'rmag': 'r-band magnitude',
    'imag': 'i-band magnitude',
    'zmag': 'z-band magnitude',
    'DDO51_finfmag': '',
    'teff': '[K] effective temperature in TRILEGAL simulation',
    'rad': '[solar radii] stellar radius',
    'Lum': '[solar luminosities] stellar luminosity',
    'GLon': '[degrees] galactic longitude',
    'GLat': '[degrees] galactic latitude',
    'Imag_reddened': '',
    'd': '',
    'ELon': '[degrees] ecliptic longitude',
    'ELat': '[degrees] ecliptic latitude',
    'tred': '',
    'region': '',
    'max_T': '',
    'Pdet_fixedBeta': 'detection probability',
    'SNR_fixedBeta': 'signal-to-noise ratio',
    'Pdet_varyBeta': 'detection probability',
    'SNR_varyBeta': 'signal-to-noise ratio',
    'P_mix': 'average detection probability',
    'Rank_Pmix': 'rank of average detection probability'
}

history_header, history_data = mesa.load_history(folder + '/LOGS/history.data')
profile_header, profile_data = mesa.load_profile(folder + '/final.profile')
gyre_header, gyre_profile = gyre.load_gyre(folder + '/final.profile.GYRE.rot')

nml = AADG3.load_namelist(basename + '.in')

# header data for FITS output
    
header = fits.Header()
header['ORIGIN'] = ('Uni. Birmingham', 'institution responsible for creating this file')
header['DATE'] = (datetime.today().strftime('%Y-%m-%d'), 'file creation date')
header['ID'] = (ID, 'ID number')
header['SECTOR'] = (sector, 'TESS observing sector')
header['SEC_RANK'] = (sec_rank, 'rank in this sector')
header['TOT_RANK'] = (int(meta['Rank_Pmix']),
                      "rank in whole sky")
header['PMIX'] = (meta['P_mix'], 'detection probability')

header['MASS'] = (history_data[-1]['star_mass'], '[solar masses] stellar mass')
header['RADIUS'] = (10.**history_data[-1]['log_R'], '[solar radii] stellar radius')
header['AGE'] = (history_data[-1]['star_age']/1e9, '[Gyr] stellar age')
header['TEFF'] = (10.**history_data[-1]['log_Teff'], '[K] effective temperature')
header['LOGG'] = (history_data[-1]['log_g'], '[cm/s2] log10 surface gravity')
header['LUM'] = (10.**history_data[-1]['log_L'], '[solar luminosities] stellar luminosity')
header['X_C'] = (profile_data[-1]['x'], 'central hydrogen abundance')
header['Y_C'] = (profile_data[-1]['y'], 'central helium abundance')
header['Z_INI'] = (history_header['initial_z'][()], 'initial metal abundance')

Zsun = 0.01756
Ysun = 0.26618
Xsun = 1.0 - Zsun - Ysun
dq = 10.**profile_data['logdq']
I = (profile_data['q'] > 0.9999)
FeH = np.log10(np.sum(profile_data['z'][I]*dq[I])/np.sum(profile_data['x'][I]*dq[I])/(Zsun/Xsun))

header['FE_H'] = (FeH, 'final metallicity [Fe/H]')

header['OMEGA_C'] = (gyre_profile['Omega'][0]/2./np.pi*1e6, '[uHz] central/core rotation rate')
header['OMEGA_E'] = (gyre_profile['Omega'][-1]/2./np.pi*1e6, '[uHz] surface/envelope rotation rate')

header['VR'] = (meta['vr'], '[km/s] radial velocity')
header['MU0'] = (meta['mu0'], 'distance modulus')
header['AV'] = (meta['Av'], 'interstellar reddening')

header['ELON'] = (meta['ELon'], '[degrees] ecliptic longitude')
header['ELAT'] = (meta['ELat'], '[degrees] ecliptic latitude')
header['GLON'] = (meta['gall'], '[degrees] galactic longitude')
header['GLAT'] = (meta['galb'], '[degrees] galactic latitude')
header['GC'] = (int(meta['Gc']), 'gal. comp.: 1-4 = thin/thick disc, halo, bulge')
header['COMP'] = (int(meta['comp']), 'binarity: 0-2 = single, primary, secondary')
header['SIGMA'] = (meta['sigma'], '[ppm] white noise amplitude')

header['SEED'] = (nml['user_seed'], 'seed for random number generator')
header['N_CADS'] = (nml['n_cadences'], 'number of cadences in hemisphere')
header['GRAN_SIG'] = (nml['sig'], '[ppm] granulation amplitude')
header['GRAN_TAU'] = (nml['tau'], '[s] granulationtimescale')
header['INC'] = (nml['inclination'], '[degrees] inclination')

# lightcurve data for FITS output

# flux = np.loadtxt(args.asc)
with open(args.asc, 'r') as f:
    flux = np.array([float(line) for line in f.readlines()])

data = np.zeros(len(flux), dtype=[('TIME', '>f8'), ('FLUX', '>f4'), ('CADENCENO', '>i4')])
data['FLUX'] = flux
data['CADENCENO'] = sector_starts[sector] + np.arange(len(flux), dtype=int)
data['TIME'] = 120.0*data['CADENCENO']/86400.0

data_fits = fits.BinTableHDU(data)
data_fits.header['TTYPE1'] = (data_fits.header['TTYPE1'], 'column title: data timestamp')
data_fits.header['TUNIT1'] = ('d', 'column unit: days')
data_fits.header['TTYPE2'] = (data_fits.header['TTYPE2'], 'column title: intensity variation')
data_fits.header['TUNIT2'] = ('ppm', 'column unit: parts per million')
data_fits.header['TTYPE3'] = (data_fits.header['TTYPE3'], 'column title: cadence number')

for i in range(1, 1 + data_fits.header['TFIELDS']):
    k = 'TFORM%i' % i
    form = data_fits.header[k]
    data_fits.header[k] = (form, forms[form])

# mode data for FITS ouptut
    
modes = np.loadtxt(basename + '.con',
                   dtype=[('L', '>i4'), ('N', '>i4'), ('FREQ', '>f8'),
                          ('WIDTH', '>f8'), ('POWER', '>f8'), ('ROT', '>f8')])

rot = np.genfromtxt(basename + '.rot', names=['n', 'l', 'm', 'dnu'])
modes['ROT'][modes['L']>0] = rot[rot['m']==1]['dnu']

modes_fits = fits.BinTableHDU(modes)
modes_fits.header['TTYPE1'] = (modes_fits.header['TTYPE1'], 'column title: angular degree')
modes_fits.header['TTYPE2'] = (modes_fits.header['TTYPE2'], 'column title: radial order')
modes_fits.header['TTYPE3'] = (modes_fits.header['TTYPE3'], 'column title: mode frequency')
modes_fits.header['TUNIT3'] = ('uHz', 'column units: microhertz')
modes_fits.header['TTYPE4'] = (modes_fits.header['TTYPE4'], 'column title: linewidth')
modes_fits.header['TUNIT4'] = ('uHz', 'column units: microhertz')
modes_fits.header['TTYPE5'] = (modes_fits.header['TTYPE5'], 'column title: mode power')
modes_fits.header['TUNIT5'] = ('ppm2', 'column units: parts per million squared')
modes_fits.header['TTYPE6'] = (modes_fits.header['TTYPE6'], 'column title: rotational splitting')
modes_fits.header['TUNIT6'] = ('uHz', 'column units: microhertz')

for i in range(1, 1 + modes_fits.header['TFIELDS']):
    k = 'TFORM%i' % i
    form = modes_fits.header[k]
    modes_fits.header[k] = (form, forms[form])

hdul = fits.HDUList([fits.PrimaryHDU(header=header), data_fits,
                     # fits.TableHDU(header=atl_fits, name='ATL'),
                     # fits.TableHDU(header=tri_fits, name='TRILEGAL'),
                     modes_fits])
hdul.writeto(fitsname, overwrite=True)

if not args.skip_test:
    test = fits.open(fitsname)
    # print(test[0].header)
    test.close()
