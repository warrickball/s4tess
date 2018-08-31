# Version 0

## Cross matching and extracting data from ATL and TRILEGAL

The first few steps only need to be carried out once after which input
for subsequent steps is always available.  These steps shoudn't
change.

    mkdir data
    python3 ../scripts/get_metadata.py ../data/atl.npy ../data/tri.npy data/meta.npy -v -N 1000

## Workflow

In this example, we do the first 20 models.

### Make folders

    mkdir alpha
    seq -w 0 00019 | xargs -I{} -t cp -r ../template/ alpha/{}

### Export metadata

    python3 ../scripts/export_metadata.py data/meta.npy "alpha/#/#.meta" -N 20 -v --fmt {:05d}

### Make MESA input

    seq -w 0 00019 | xargs -t -I{} python3 ../scripts/make_MESA_input.py alpha/{}/{}.meta --Tc ../data/Tc.dat

or

    ls alpha/*/*.meta | xargs -t -I{} python3 ../scripts/make_MESA_input.py {} --Tc ../data/Tc.dat

### Run MESA

    sbatch --export=RUN="alpha" --array=0-19 MESA_sector_array.slurm

Skips folders in which ``final.profile.GYRE`` exists.

### Add rotation

    seq -w 0 00019 | xargs -t -I{} bash add_rotation.sh alpha/{}

### Run GYRE

    sbatch --export=RUN="alpha" --array=0-19 GYRE_sector_array.slurm

Skips folders in which ``gyre_summary.txt`` exists.

### Make AADG3 input

    seq -w 0 00019 | xargs -t -I{} python3 ../scripts/make_AADG3_input.py alpha/{}

### Run AADG3

    sbatch --export=RUN="alpha" --array=0-19 AADG3_sector_array.slurm

Skips folders in which ``"$ID".asc`` exists.

### Add white noise

    seq -w 0 00019 | xargs -t -I{} python3 ../scripts/add_white_noise.py alpha/{}

### Separate data into sectors

    seq -w 0 00019 | xargs -t -I{} python3 ../scripts/sectorize.py alpha/{}

Skips folders where ``"$ID"_WN.asc`` doesn't exist.

### Convert output to FITS

    ls alpha/000[01][0-9]/*_WN_*.asc | xargs -t -n1 python3 ../scripts/asc_to_fits.py

### Save power spectra

    ls alpha/000*/*_WN_*.asc | xargs -n1 dirname | xargs -t -n1 bash save_AADG3_PS.sh

