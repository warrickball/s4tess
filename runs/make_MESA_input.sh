#!/usr/bin/env bash

# Does all calculations up to prepare data up to running MESA.

if [ -z "$3" ]
then
    echo "Usage: bash make_MESA_input.sh <hemisphere> <start star> <end star>"
    echo "       where <hemisphere> is either 'north' or 'south'"
    exit 1
fi

TOP=..

for STAR in $(seq -w $2 $3)
do
    mkdir -p $1/$STAR/LOGS
    mkdir -p $1/$STAR/photos
    cp -n $TOP/template/* "$1"/$STAR/
done
python3 $TOP/scripts/make_MESA_input.py data/tri_"$1".npy "$1"/{:05d} -v --Tc $TOP/data/Tc.dat -N $2 $3
python3 $TOP/scripts/extract_catalogue_data.py data/atl_"$1".npy "$1"/{:05d}/data.atl -N $2 $3
python3 $TOP/scripts/extract_catalogue_data.py data/tri_"$1".npy "$1"/{:05d}/data.tri -N $2 $3
for STAR in $(seq -w $2 $3)
do
    ln -s $(cd $TOP ; pwd)/mesa_work/star "$1"/$STAR/star
    mv "$1"/$STAR/data.atl "$1"/$STAR/"$STAR".atl
    mv "$1"/$STAR/data.tri "$1"/$STAR/"$STAR".tri
done

echo "Done."
