#!/usr/bin/env bash

# given a folder, create AADG3 input for one run

OLDPWD=$(pwd)
# BASENAME=$(echo $1 | sed 's:/:_:g')
BASENAME=$(basename $1)
cd $1
python3 ../../../../scripts/make_AADG3_input.py gyre_summary.txt $BASENAME.in $BASENAME.con $BASENAME.rot --min-H 0.0001
cd $OLDPWD
