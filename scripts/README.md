# Scripts

## Workflow

### 1. Gather TRILEGAL output

    python3 combine_trilegal.py ../data/tri.npy

Combines the files in default locations into `../data/tri.npy`.

### 2. Get best targets in each sector from ATL

    python3 get_sector_best.py ../data/atl.npy ../data/atl_best1000_{:2d}.npy -N 1000

Goes down the likelihood-sorted ATL file `atl.npy` recording in which
sector each star occurs, then writes the best targets found in this
way to files with names patching `atl_best1000_{:2d}.npy`.

### 3. Find details of best targets in TRILEGAL data

    python3 cross_match_trilegal_atl.py ../data/tri.npy ../data/atl_best1000_13.npy -o ../data/tri_best1000_13.npy

Having found the best targets using the ATL data, we now go back and
find out which stars these were in the original TRILEGAL data.  The
example would give sector 13 (sector 0 in Southern hemisphere).

### 4. Create MESA input for each star in each sector

### 5. ???

### 6. Profit
