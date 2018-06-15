#!/usr/bin/env bash

# given a folder, adds rotation to the final model if a rotating model
# doesn't already exist

if [ ! -f $1/final.profile.GYRE ]
then
    echo "ERROR: $1/final.profile.GYRE doesn't exist!"
    exit 1
fi

OLDPWD=$(pwd)
if [ ! -f $1/final.profile.GYRE.rot ]
then
    cd $1
    python3 ../../../../scripts/add_rotation.py final.profile.GYRE LOGS/history.data --no-env-noise
    cd $OLDPWD
fi
