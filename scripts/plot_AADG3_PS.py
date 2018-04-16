#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
import AADG3
from astropy.stats import LombScargle as LS
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('AADG3_input', type=str,
                    help="filename of AADG3 namelist")
args = parser.parse_args()

nml, con, rot = AADG3.load_all_input(args.AADG3_input)
y = np.loadtxt(nml['nameout'])
t = np.arange(len(y), dtype=float)*nml['cadence']
f, p = LS(t, y).autopower(normalization='psd', nyquist_factor=1.0,
                          samples_per_peak=1)
f = f*1e6  # to uHz
p = p*np.mean(y**2)/np.sum(p)/(f[1]-f[0])  # Bill's normalization

pl.loglog(f, p, label='time series')
pl.loglog(f, AADG3.PS_model(f, nml, con, rot), label='model')
pl.legend()
pl.show()
