#!/usr/bin/env python

import numpy as np
from tomso import mesa
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
args = parser.parse_args(['../runs/v0/south/00000/00000_WN_11_0000.asc'])

forms = {'J': 'column format: signed 32-bit integer',
         'D': 'column format: 64-bit floating point',
         'E': 'column format: 32-bit floating point'}
cads_per_sector = 720*137//5  # (cadences/day)*(days/sector)

ascname = args.asc.split('/')[-1]
folder = '/'.join(args.asc.split('/')[:-1])
basename = folder + '/' + args.asc.split('/')[-2]
sector = int(ascname.split('.')[0].split('_')[2])
# if 'north' in args.asc:
#     sector += 13

if args.output:
    fitsname = args.output
else:
    if args.asc.endswith('.asc'):
        fitsname = args.asc[:-4] + '.fits'
    else:
        fitsname = args.asc + '.fits'

tri = {}  # TRILEGAL data
with open(basename + '.tri', 'r') as f:
    for line in f.readlines():
        k, v = line.split('=')
        tri[k.strip()] = float(v)

tri_comments = {
    'gall': 'galactic longitude',
    'galb': 'galactic latitude',
    'Gc': 'galactic component: 1--4 = thin disc, thick disc, halo, bulge',
    'logAge': '[yr] log10 age',
    'M_H': '',
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
    'TESSmag': '',
    'Jmag': '',
    'Hmag': '',
    'Ksmag': '',
    'Keplermag': '',
    'gmag': '',
    'rmag': '',
    'imag': '',
    'zmag': '',
    'DDO51_finfmag': ''
}

atl = {}  # ATL data
with open(basename + '.atl', 'r') as f:
    for line in f.readlines():
        k, v = line.split('=')
        atl[k.strip()] = float(v)

atl_comments = {
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
    'numax': '[uHz] frequency of maximum oscillation power',
    'max_T': '',
    'Pdet_fixedBeta': 'detection probability',
    'SNR_fixedBeta': 'signal-to-noise ratio',
    'Pdet_varyBeta': 'detection probability',
    'SNR_varyBeta': 'signal-to-noise ratio',
    'P_mix': 'average detection probability',
    'Rank_Pmix': 'rank of average detection probability'
}

v = np.loadtxt(basename + '.vr')[()]

history_header, history_data = mesa.load_history(folder + '/LOGS/history.data')
profile_header, profile_data = mesa.load_profile(folder + '/final.profile.GYRE')

# header data for FITS output
    
header = fits.Header()
header['ORIGIN'] = ('Uni. Birmingham', 'institution responsible for creating this file')
header['DATE'] = (datetime.today().strftime('%Y-%m-%d'), 'file creation date')
header['TEFF'] = (10.**history_data[-1]['log_Teff'], '[K] Effective temperature')
header['LOGG'] = (history_data[-1]['log_g'], '[cm/s2] log10 surface gravity')
header['SECTOR'] = (sector, 'Observing sector')
header['RADIUS'] = (10.**history_data[-1]['log_R'], '[solar radii] stellar radius')
header['V_R'] = (v, '[km/s] radial velocity')

# lightcurve data for FITS output

flux = np.loadtxt(args.asc)
data = np.zeros(len(flux), dtype=[('TIME', '>f8'), ('FLUX', '>f4'), ('CADENCENO', '>i4')])
data['FLUX'] = flux
data['CADENCENO'] = sector*cads_per_sector + np.arange(len(flux), dtype=int)
data['TIME'] = 120.0*data['CADENCENO']/86400.0

data_fits = fits.BinTableHDU(data)
data_fits.header['TTYPE1'] = (data_fits.header['TTYPE1'], 'column title: data timestamp')
data_fits.header['TUNIT1'] = ('d', 'column unit: days')
data_fits.header['TTYPE2'] = (data_fits.header['TTYPE2'], 'column title: intensity variation')
data_fits.header['TUNIT2'] = ('', '')
data_fits.header['TTYPE3'] = (data_fits.header['TTYPE3'], 'column title: cadence number')

for i in range(1, 1 + data_fits.header['TFIELDS']):
    k = 'TFORM%i' % i
    form = data_fits.header[k]
    data_fits.header[k] = (form, forms[form])

# ATL data for FITS ouptut

atl_fits = {k: (atl[k], atl_comments[k]) for k in atl.keys()}
atl_fits = fits.Header()
for k in atl.keys():
    atl_fits[k] = (atl[k], atl_comments[k])
    
# TRILEGAL data for FITS ouptut
    
tri_fits = {k: (tri[k], tri_comments[k]) for k in tri.keys()}
tri_fits = fits.Header()
for k in tri.keys():
    tri_fits[k] = (tri[k], tri_comments[k])

# mode data for FITS ouptut
    
modes = np.loadtxt(basename + '.con',
                   dtype=[('l', '>i4'), ('n', '>i4'), ('nu', '>f8'),
                          ('width', '>f8'), ('amp2', '>f8'), ('rot', '>f8')])

rot = np.genfromtxt(basename + '.rot', names=['n', 'l', 'm', 'dnu'])
modes['rot'][modes['l']>0] = rot[rot['m']==1]['dnu']

modes_fits = fits.BinTableHDU(modes)
modes_fits.header['TTYPE1'] = (modes_fits.header['TTYPE1'], 'column title: angular degree')
modes_fits.header['TTYPE2'] = (modes_fits.header['TTYPE2'], 'column title: radial order')
modes_fits.header['TTYPE3'] = (modes_fits.header['TTYPE3'], 'column title: mode frequency')
modes_fits.header['TUNIT3'] = ('uHz', 'column units: microhertz')
modes_fits.header['TTYPE4'] = (modes_fits.header['TTYPE4'], 'column title: linewidth')
modes_fits.header['TUNIT4'] = ('uHz', 'column units: microhertz')
modes_fits.header['TTYPE5'] = (modes_fits.header['TTYPE5'], 'column title: mode power')
modes_fits.header['TTYPE6'] = (modes_fits.header['TTYPE6'], 'column title: rotational splitting')
modes_fits.header['TUNIT6'] = ('uHz', 'column units: microhertz')

for i in range(1, 1 + modes_fits.header['TFIELDS']):
    k = 'TFORM%i' % i
    form = modes_fits.header[k]
    modes_fits.header[k] = (form, forms[form])

hdul = fits.HDUList([fits.PrimaryHDU(header=header), data_fits,
                     fits.TableHDU(header=atl_fits, name='ATL'),
                     fits.TableHDU(header=tri_fits, name='TRILEGAL'),
                     modes_fits])
hdul.writeto(fitsname, overwrite=True)

test = fits.open(fitsname)
print(test[0].header)
# test.close()
