# Scripts

## Workflow

### 1. Gather TRILEGAL output

    python3 combine_trilegal.py ../data/tri.npy

Combines the files in default locations into `../data/tri.npy`.  Much
easier and faster to manipulate one NumPy binary file.

### 2. Convert TRILEGAL data to CSV data for Mat to run ATL

    python3 trilegal_to_atl.py ../data/tri.npy ~/mnt/adfbison/data/trilegal_for_atl.dat

Extracts critical data from TRILEGAL simulation that Mat can use as
input for ATL.

### 3. Convert ATL CSV output to NumPy binary

    python3 atl_csv_to_npy.py ~/mnt/adfbison/data/trilegal_results.csv ../data/atl.npy

### 4. Get best targets in each sector from ATL output

    python3 get_sector_best.py ../data/atl.npy ../data/atl_best1000_{:2d}.npy -N 1000

Goes down the likelihood-sorted ATL file `atl.npy` recording in which
sector each star occurs, then writes the best targets found in this
way to files with names patching `atl_best1000_{:2d}.npy`.

### 5. Find details of best targets in TRILEGAL data

    python3 cross_match_trilegal_atl.py ../data/tri.npy ../data/atl_best1000_13.npy -o ../data/tri_best1000_13.npy

Having found the best targets using the ATL data, we now go back and
find out which stars these were in the original TRILEGAL data.  The
example would give sector 13 (sector 0 in Southern hemisphere).

### 6. Create MESA input for each star in each sector

    python3 make_MESA_input.py ../data/tri_best1000_13.npy ../results/13/{:04d}/inlist_run

Loops through all the data in first argument and writes to second,
with star number inserted.

### 7. ???

### 8. Profit
