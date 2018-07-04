#!/usr/bin/env bash

# given a timeseries, save power spectrum, with AADG3 model

OLDPWD=$(pwd)
# BASENAME=$(echo $1 | sed 's:/:_:g')
DIRNAME=$(dirname $1)
FILENAME=$(basename $1 .asc)
ID=$(basename $DIRNAME)
cd $DIRNAME
python3 ../../../../scripts/save_AADG3_PS.py "$ID".in "$FILENAME".pow --nameout "$FILENAME".asc
cd $OLDPWD
