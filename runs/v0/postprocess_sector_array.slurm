#!/usr/bin/env bash

# when using this command, you must export the SECTOR variable
# and specify the array range
# e.g. sbatch --export=SECTOR="north" --array=0-99 AADG3_sector_array.slurm

#SBATCH --nodes 1
#SBATCH --tasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 256M
#SBATCH --time 00:10:00
#SBATCH --output slurm/postprocess_sector_array_%A_%a.out
#SBATCH --error slurm/postprocess_sector_array_%A_%a.err


module purge; module load bluebear  # this line is required
module load scipy/0.19.0-python-3.5.2

echo $SECTOR
ID=$(printf "%05d" $SLURM_ARRAY_TASK_ID)
cd $SLURM_SUBMIT_DIR
python3 ../../scripts/add_white_noise.py $SECTOR/$ID
python3 ../../scripts/sectorize.py $SECTOR/$ID
ls $SECTOR/$ID/"$ID"_WN_*.asc | xargs -t -n1 python3 ../../scripts/asc_to_fits.py
