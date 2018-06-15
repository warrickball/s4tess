#!/usr/bin/env bash

# remember: MESA takes about 15 core minutes per star
# when using this command, you must export the SECTOR variable
# e.g. sbatch --export=SECTOR=00 mesa_sector_array.slurm

#SBATCH --nodes 1
#SBATCH --tasks 1
#SBATCH --mem 1G
#SBATCH --time 00:10:00
#SBATCH --output slurm/AADG3_sector_array_%A_%a.out
#SBATCH --error slurm/AADG3_sector_array_%A_%a.err
#SBATCH --array 0-399

module purge; module load bluebear  # this line is required

export OMP_NUM_THREADS=$SLURM_NTASKS
echo OMP_NUM_THREADS $OMP_NUM_THREADS
echo $SECTOR
ID=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
cd $SLURM_SUBMIT_DIR/$SECTOR/$ID
pwd
if [ -f "$SECTOR"_"$ID".asc ]
then
    echo "ERROR: found $(pwd)/${SECTOR}_${ID}.asc"
    echo "       AADG3 appears to have already finished successfully."
    exit 1
fi

AADG3 "$SECTOR"_"$ID".in
