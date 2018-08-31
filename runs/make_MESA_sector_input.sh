#!/usr/bin/env bash

# Does all calculations up to prepare data up to running MESA.

if [ -z "$3" ]
then
    echo "Usage: bash make_MESA_input.sh <start sector> <end sector> <number of stars>"
    exit 1
fi

NSECTORS=$(seq $1 $2 | wc -l)
NSTARS=$3
TOP=..

for SECTOR in $(seq -w $1 $2)
do
    echo "Preparing sector $SECTOR of $NSECTORS..."
    mkdir -p $SECTOR
    for STAR in $(seq -w 0000 $((NSTARS-1)))
    do
	mkdir -p $SECTOR/$STAR
	cp $TOP/template/* $SECTOR/$STAR/
    done
    python3 $TOP/scripts/make_MESA_input.py data/tri_v0_$SECTOR.npy $SECTOR/{:04d} -v --Tc $TOP/data/Tc.dat -N $NSTARS
    python3 $TOP/scripts/extract_catalogue_data.py data/atl_v0_$SECTOR.npy $SECTOR/{:04d}/data.atl -N $NSTARS
    python3 $TOP/scripts/extract_catalogue_data.py data/tri_v0_$SECTOR.npy $SECTOR/{:04d}/data.tri -N $NSTARS
    for STAR in $(seq -w 0000 $((NSTARS-1)))
    do
	ln -s $(cd $TOP ; pwd)/mesa_work/star $SECTOR/$STAR/star
	mv $SECTOR/$STAR/data.atl $SECTOR/$STAR/"$SECTOR"_"$STAR".atl
	mv $SECTOR/$STAR/data.tri $SECTOR/$STAR/"$SECTOR"_"$STAR".tri
    done
done

echo "Done."
