# Version 0

## Cross matching and extracting data from ATL and TRILEGAL

The first few steps only need to be carried out once after which input
for subsequent steps is always available.  These steps shoudn't
change.

    mkdir data
    python3 ../../scripts/get_sector_best.py ../../data/atl.npy data/atl_v0_{:02d}.npy -N 1000 -v
    seq -w 00 25 | xargs -t -I{} python3 ../../scripts/cross_match_trilegal_atl.py ../../data/tri.npy data/atl_v0_{}.npy -o data/tri_v0_{}.npy -v

## Workflow

### Make main folders with MESA input

    bash make_MESA_input.sh

Cross matches catalogues and generates priority target lists.

### Run MESA

    seq -w 00 25 | xargs -I{} sbatch --export=SECTOR={} mesa_sector_array.slurm

Skips folders in which ``final.profile.GYRE`` exists.

### Add rotation

    ls */*/final.profile.GYRE | xargs -n1 dirname | xargs -t -n1 bash add_rotation.sh

### Run GYRE

    seq -w 00 25 | xargs -I{} sbatch --export=SECTOR={} gyre_sector_array.slurm

Skips folders in which ``gyre_summary.txt`` exists.

### Make AADG3 input

    ls */*/gyre_summary.txt | xargs -n1 dirname | xargs -t -n1 bash make_AADG3_input.sh

### Run AADG3

    seq -w 00 25 | xargs -I{} sbatch --export=SECTOR={} AADG3_sector_array.slurm

Skips folders in which ``"$SECTOR"_"$ID".asc`` exists.

### Add white noise

    ls */*/*.asc | xargs -n1 dirname | xargs -t -n1 bash add_white_noise.sh

### Save power spectra

    ls */*/*_WN.asc | xargs -n1 dirname | xargs -t -n1 bash save_AADG3_PS.sh
