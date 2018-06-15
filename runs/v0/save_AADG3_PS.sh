#!/usr/bin/env bash

# given a folder, create AADG3 input for one run

OLDPWD=$(pwd)
BASENAME=$(echo $1 | sed 's:/:_:g')
cd $1
python3 ../../../../scripts/save_AADG3_PS.py $BASENAME.in $BASENAME.pow
cd $OLDPWD
