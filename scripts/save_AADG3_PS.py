#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
import AADG3
from astropy.stats import LombScargle
from argparse import ArgumentParser

def LS(t, y):
    f, p = LombScargle(t, y).autopower(normalization='psd',
                                       nyquist_factor=1.0, samples_per_peak=1)
    f = f*1e6  # to uHz
    p = p*np.mean(y**2)/np.sum(p)/(f[1]-f[0])  # Bill's normalization
    return f, p

parser = ArgumentParser()
parser.add_argument('AADG3_input', type=str,
                    help="filename of AADG3 namelist")
parser.add_argument('output', type=str, help="output file")
parser.add_argument('--fmt', type=str, default='%16.7e',
                    help="format statement for numpy.savetxt")
args = parser.parse_args()

nml, con, rot = AADG3.load_all_input(args.AADG3_input)
try:
    filename = nml['nameout'].replace('.asc', '_WN.asc')
    y = np.loadtxt(filename)
except OSError:
    filename = '/'.join(args.AADG3_input.split('/')[:-1]) \
               + '/' + nml['nameout'].replace('.asc', '_WN.asc')
    y = np.loadtxt(filename)
    
t = np.arange(len(y), dtype=float)*nml['cadence']
f, p = LS(t, y)
p0 = AADG3.PS_model(f, nml, con, rot)

np.savetxt(args.output, np.vstack((f, p, p0)).T, fmt=args.fmt)
