# Version 0

## Cross matching and extracting data from ATL and TRILEGAL

The first few steps only need to be carried out once after which input
for subsequent steps is always available.  These steps shoudn't
change.

    mkdir data
    python3 ../../scripts/get_sector_best.py ../../data/atl.npy data/atl_v0_{:02d}.npy -N 1000 -v
    seq -w 00 25 | xargs -t -I{} python3 ../../scripts/cross_match_trilegal_atl.py ../../data/tri.npy data/atl_v0_{}.npy -o data/tri_v0_{}.npy -v

## Get unique targets in Northern and Southern hemispheres

    python3 unique.py north
    python3 unique.py south

## Workflow

In this example, we do the first 100 models in the Northern hemisphere.

### Make main folders with MESA input

    bash make_MESA_input.sh north 00000 00099

### Run MESA

    sbatch --export=SECTOR="north" --array=0-99 MESA_sector_array.slurm

Skips folders in which ``final.profile.GYRE`` exists.

### Add rotation

    ls north/*/final.profile.GYRE | xargs -n1 dirname | xargs -t -n1 bash add_rotation.sh

### Run GYRE

    sbatch --export=SECTOR="north" --array=0-99 GYRE_sector_array.slurm

Skips folders in which ``gyre_summary.txt`` exists.

### Make AADG3 input

    ls north/*/gyre_summary.txt | xargs -n1 dirname | xargs -t -n1 python3 ../../scripts/make_AADG3_input.py

### Run AADG3

    sbatch --export=SECTOR="north" --array=0-99 AADG3_sector_array.slurm

Skips folders in which ``"$ID".asc`` exists.

### Add white noise

    ls north/*/*.asc | xargs -n1 dirname | xargs -t -n1 bash add_white_noise.sh

### Separate data into sectors

    python3 sectorize.py north

Skips folders where ``"$ID"_WN.asc`` doesn't exist.

### Convert output to FITS

    ls north/*/*_WN_*.asc | xargs -t -n1 python3 ../../scripts/asc_to_fits.py

### Save power spectra

    ls north/*/*_WN_*.asc | xargs -n1 dirname | xargs -t -n1 bash save_AADG3_PS.sh

