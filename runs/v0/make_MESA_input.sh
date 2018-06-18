#!/usr/bin/env bash

# Does all calculations up to prepare data up to running MESA.

if [ -z "$2" ]
then
    echo "Usage: bash make_MESA_input.sh <start star> <end star>"
    exit 1
fi

TOP=../..

for STAR in $(seq -w $1 $2)
do
    mkdir -p unique/$STAR
    cp -n $TOP/template/* unique/$STAR/
done
python3 $TOP/scripts/make_MESA_input.py data/tri_unique.npy unique/{:05d} -v --Tc $TOP/data/Tc.dat -N $1 $2
python3 $TOP/scripts/extract_catalogue_data.py data/atl_unique.npy unique/{:05d}/data.atl -N $1 $2
python3 $TOP/scripts/extract_catalogue_data.py data/tri_unique.npy unique/{:05d}/data.tri -N $1 $2
for STAR in $(seq -w $1 $2)
do
    ln -s $(cd $TOP ; pwd)/mesa_work/star unique/$STAR/star
    mv unique/$STAR/data.atl unique/$STAR/"$STAR".atl
    mv unique/$STAR/data.tri unique/$STAR/"$STAR".tri
done

echo "Done."
