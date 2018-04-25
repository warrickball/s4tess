#!/usr/bin/env python

import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser(description=
"""Add white noise to a lightcurve using the formula in Stassun et
al. (2017) for an I-band magnitude in the given file.""")
parser.add_argument('lc_in', type=str,
                    help='filename of input lightcurve')
parser.add_argument('Imag', type=str,
                    help="I-band magnitude or text file with a line"
                    "containing 'imag = ...'")
parser.add_argument('lc_out', type=str,
                    help='filename of output lightcurve')
args = parser.parse_args()

try:
    Imag = float(args.Imag)
except ValueError:
    with open(args.Imag, 'r') as f:
        lines = [line.strip().split('=') for line in f.readlines()]
        for line in lines:
            if line[0].lower().startswith('imag'):
                Imag = float(line[1])
                break
            else:
                raise IOError("no line imag = ... in file %s" % args.Imag)
        

# see Stassun et al. (2017)
# https://arxiv.org/abs/1706.00495

# if Imag is less than local minimum aroung 5,
# use flat line that roughly corresponds to systematic noise
Imag = np.maximum(Imag, 4.97756051)
F = 4.73508403525e-5
E = -0.0022308015894
D = 0.0395908321369
C = -0.285041632435
B = 0.850021465753
lnA = 3.29685004771
lnsigma = lnA + B*Imag + C*Imag**2 + D*Imag**3 + E*Imag**4 + F*Imag**5
sigma = np.exp(lnsigma)
        
asc = np.loadtxt(args.lc_in)
asc = asc + np.random.randn(len(asc))*sigma
np.savetxt(args.lc_out, asc)
