#!/usr/bin/env bash

# bash script to prepare data for a full calculation of the mock
# catalogue, starting with the NumPy files for the ATL's top-ranked
# targets

TOP_DIR=../..

export OMP_NUM_THREADS=4
CORES=4
RUNS=$(echo $CORES/$OMP_NUM_THREADS | bc)
NSTARS=1
NSECTORS=1

function do_one() {
    ONE_TOP_DIR=../../../..
    cd $1/$2
    ln -s $ONE_TOP_DIR/mesa_work/star
    MAX_TRIES=10
    TRIES=1
    echo TRY $TRIES > mesa.log
    ./rn >> mesa.log
    while grep --quiet convergence mesa.log && ((TRIES <= MAX_TRIES))
    do
	python3 $ONE_TOP_DIR/scripts/modify_Tc.py inlist_run 5000 -d
	((TRIES++))
	echo TRY $TRIES > mesa.log
	./rn >> mesa.log
    done
    python3 $ONE_TOP_DIR/scripts/add_rotation.py final.profile.GYRE LOGS/history.data
    $GYRE_DIR/bin/gyre gyre.in > gyre.log
    python3 $ONE_TOP_DIR/scripts/make_AADG3_input.py gyre_summary.txt $1_$2.in $1_$2.con $1_$2.rot
    AADG3 $1_$2.in
    python3 $ONE_TOP_DIR/scripts/add_white_noise.py atl_data.txt tri_data.txt $1_$2.asc $1_$2_WN.asc
}

echo
date
echo

echo "OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo "CORES           = $CORES"
echo "RUNS            = $RUNS"
echo "NSTARS          = $NSTARS"
echo "NSECTORS        = $NSECTORS"
echo

export -f do_one

seq -w 00 99 | xargs rm -rf
test -d data || mkdir data
rm -f data/*_test_*.npy

python3 $TOP_DIR/scripts/get_sector_best.py $TOP_DIR/data/atl.npy data/atl_test_{:02d}.npy -N $NSTARS -v
for SECTOR in $(seq -w 00 $((NSECTORS-1)))
do
    echo
    echo "-----------"
    echo " SECTOR $SECTOR"
    echo "-----------"
    echo

    python3 $TOP_DIR/scripts/cross_match_trilegal_atl.py $TOP_DIR/data/tri.npy data/atl_test_$SECTOR.npy -o data/tri_test_$SECTOR.npy -v
    mkdir -p $SECTOR
    seq -w 0000 $((NSTARS-1)) | xargs -t -I{} cp -R $TOP_DIR/template $SECTOR/{}
    python3 $TOP_DIR/scripts/make_MESA_input.py data/tri_test_$SECTOR.npy $SECTOR/{:04d} -v --Tc $TOP_DIR/data/Tc.dat
    python3 $TOP_DIR/scripts/extract_catalogue_data.py data/atl_test_$SECTOR.npy $SECTOR/{:04d}/atl_data.txt
    python3 $TOP_DIR/scripts/extract_catalogue_data.py data/tri_test_$SECTOR.npy $SECTOR/{:04d}/tri_data.txt
    seq -w 0000 $((NSTARS-1)) | xargs -t -P$RUNS -I{} bash -c "do_one $SECTOR {}"
done

echo
date
echo
