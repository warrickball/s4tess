#!/usr/bin/env python3

import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('csv', type=str,
                    help="input CSV file with ATL results")
parser.add_argument('npy', type=str,
                    help="output NumPy binary with ATL results")
args = parser.parse_args()

data = np.genfromtxt(args.csv, delimiter=',', names=True)
np.save(args.npy, data)
