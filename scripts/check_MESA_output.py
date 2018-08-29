#!/usr/bin/env python3

import numpy as np
from tomso import mesa
from argparse import ArgumentParser

def one_line(name, x, y, expx=False, expy=False):
    xx = 10.**x if expx else x
    yy = 10.**y if expy else y
    print('{:>12s}{:>12.5f}{:>12.5f}{:>12.5f}{:>12.5f}'.format(
        name, xx, yy, xx-yy, xx/yy-1.))

parser = ArgumentParser(description=
"""Cross check MESA output against data in TRILEGAL output.""")
parser.add_argument('folders', type=str, nargs='+',
                    help="folders containing results")
args = parser.parse_args()

for folder in args.folders:
    basename = folder.split('/')[-1]

    header, history = mesa.load_history(folder + '/LOGS/history.data')

    with open('%s/%s.meta' % (folder, basename), 'r') as f:
        lines = [line.strip().split(' = ') for line in f.readlines()]
        meta = {line[0]: float(line[1]) for line in lines}

    columns = ['#           ', 'TRILEGAL', 'MESA', 'diff', 'frac diff']
    header_fmt = '{:>12s}'*len(columns)
    print(header_fmt.format(*columns))

    one_line('Teff/K', meta['logTe'], history[-1]['log_Teff'], expx=True, expy=True)
    one_line('log(L/Lsun)', meta['logL'], history[-1]['log_L'])
    one_line('age/Gyr', meta['logAge']-9.0, history[-1]['star_age']/1e9, expx=True)

    print()
