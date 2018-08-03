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

atl = {}  # ATL data
with open(basename + '.atl', 'r') as f:
    for line in f.readlines():
        k, v = line.split('=')
        atl[k.strip()] = float(v)

v = np.loadtxt(basename + '.vr')[()]

history_header, history_data = mesa.load_history(folder + '/LOGS/history.data')
profile_header, profile_data = mesa.load_profile(folder + '/final.profile.GYRE')
    
header = fits.Header()
header['ORIGIN'] = ('Uni. Birmingham', 'institution responsible for creating this file')
header['DATE'] = (datetime.today().strftime('%Y-%m-%d'), 'file creation date')
header['TEFF'] = (10.**history_data[-1]['log_Teff'], '[K] Effective temperature')
header['LOGG'] = (history_data[-1]['log_g'], '[cm/s2] log10 surface gravity')
header['SECTOR'] = (sector, 'Observing sector')
header['RADIUS'] = (10.**history_data[-1]['log_R'], '[solar radii] stellar radius')
header['V_R'] = (v, '[km/s] radial velocity')

flux = np.loadtxt(args.asc)
data = np.zeros(len(flux), dtype=[('TIME', '>f8'), ('FLUX', '>f4'), ('CADENCENO', '>i4')])
data['FLUX'] = flux
data['CADENCENO'] = sector*19440 + np.arange(len(flux), dtype=int)
data['TIME'] = 120.0*data['CADENCENO']/86400.0

modes = np.loadtxt(basename + '.con',
                   dtype=[('l', '>i4'), ('n', '>i4'), ('nu', '>f8'),
                          ('width', '>f8'), ('amp2', '>f8'), ('rot', '>f8')])

rot = np.genfromtxt(basename + '.rot', names=['n', 'l', 'm', 'dnu'])
modes['rot'][modes['l']>0] = rot[rot['m']==1]['dnu']

hdul = fits.HDUList([fits.PrimaryHDU(header=header), fits.BinTableHDU(data),
                     fits.TableHDU(header=fits.Header(atl), name='ATL'),
                     fits.TableHDU(header=fits.Header(tri), name='TRILEGAL'),
                     fits.BinTableHDU(modes)])
hdul.writeto(fitsname, overwrite=True)

test = fits.open(fitsname)
print(test[0].header)
test.close()
