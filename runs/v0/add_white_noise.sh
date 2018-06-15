#!/usr/bin/env bash

# given a folder, create AADG3 input for one run

OLDPWD=$(pwd)
BASENAME=$(echo $1 | sed 's:/:_:g')
cd $1
python3 ../../../../scripts/add_white_noise.py $BASENAME.atl $BASENAME.tri $BASENAME.asc "$BASENAME"_WN.asc
cd $OLDPWD
