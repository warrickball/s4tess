#!/usr/bin/env bash

# bash script to prepare data for a full calculation of the mock
# catalogue, starting with the NumPy files for the ATL's top-ranked
# targets

export OMP_NUM_THREADS=4
CORES=8
RUNS=$(echo $CORES/$OMP_NUM_THREADS | bc)
NSTARS=100
NSECTORS=1
MAXTRIES=10

function do_one() {
    cd test_run/$1/$2
    ln -s ../../../mesa_work/star
    ./rn > mesa.log
    TRIES=1
    while grep --quiet convergence mesa.log && ((TRIES <= MAXTRIES))
    do
	python3 ../../../scripts/modify_Tc.py inlist_run 5000 -d
	./rn > mesa.log
	((TRIES++))
    done
    python3 ../../../scripts/add_rotation.py final.profile.GYRE LOGS/history.data
    $GYRE_DIR/bin/gyre gyre.in > gyre.log
    python3 ../../../scripts/make_AADG3_input.py gyre_summary.txt $1_$2.in $1_$2.con $1_$2.rot
    AADG3 $1_$2.in
    python3 ../../../scripts/add_white_noise.py atl_data.txt tri_data.txt $1_$2.asc $1_$2_WN.asc
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

rm -rf test_run/
rm data/*_test_*.npy

python3 scripts/get_sector_best.py data/atl.npy data/atl_test_{:02d}.npy -N $NSTARS -v
for sector in $(seq -w 00 $(echo $NSECTORS-1 | bc))
do
    echo
    echo "-----------"
    echo " SECTOR $sector"
    echo "-----------"
    echo

    python3 scripts/cross_match_trilegal_atl.py data/tri.npy data/atl_test_$sector.npy -o data/tri_test_$sector.npy -v
    mkdir -p test_run/$sector
    seq -w 0000 $(echo $NSTARS-1 | bc) | xargs -t -I{} cp -R template test_run/$sector/{}
    python3 scripts/make_MESA_input.py data/tri_test_$sector.npy test_run/$sector/{:04d} -v --Tc data/Tc.dat
    python3 scripts/extract_catalogue_data.py data/atl_test_$sector.npy test_run/$sector/{:04d}/atl_data.txt
    python3 scripts/extract_catalogue_data.py data/tri_test_$sector.npy test_run/$sector/{:04d}/tri_data.txt
    seq -w 0000 $(echo $NSTARS-1 | bc) | xargs -t -P$RUNS -I{} bash -c "do_one $sector {}"
done

echo
date
echo
