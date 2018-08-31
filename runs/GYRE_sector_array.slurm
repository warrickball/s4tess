#!/usr/bin/env bash

# when using this command, you must export the SECTOR variable
# and specify the array range
# e.g. sbatch --export=SECTOR="north" --array=0-99 GYRE_sector_array.slurm

#SBATCH --nodes 1
#SBATCH --tasks 1
#SBATCH --mem 1G
#SBATCH --time 00:30:00
#SBATCH --output slurm/GYRE_sector_array_%A_%a.out
#SBATCH --error slurm/GYRE_sector_array_%A_%a.err

module purge; module load bluebear  # this line is required

export OMP_NUM_THREADS=1
echo OMP_NUM_THREADS $OMP_NUM_THREADS
echo $SECTOR
export GYRE_DIR=/rds/homes/b/ballwh/code/gyre/fork
ID=$(printf "%05d" $SLURM_ARRAY_TASK_ID)
cd $SLURM_SUBMIT_DIR/$SECTOR/$ID
pwd
if [ -f gyre_summary.txt ]
then
    echo "ERROR: found $(pwd)/gyre_summary.txt"
    echo "       GYRE appears to have already finished successfully."
    exit 1
fi

$GYRE_DIR/bin/gyre gyre.in > gyre.log