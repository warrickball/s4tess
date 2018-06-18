#!/usr/bin/env bash

# remember: MESA takes about 15 core minutes per star
# when using this command, you must export the SECTOR variable
# e.g. sbatch --export=SECTOR=00 mesa_sector_array.slurm

#SBATCH --nodes 1
#SBATCH --tasks 1
#SBATCH --mem 1G
#SBATCH --time 01:30:00
#SBATCH --output slurm/AADG3_sector_array_%A_%a.out
#SBATCH --error slurm/AADG3_sector_array_%A_%a.err
#SBATCH --array 0-199

module purge; module load bluebear  # this line is required

echo $SECTOR
ID=$(printf "%05d" $SLURM_ARRAY_TASK_ID)
cd $SLURM_SUBMIT_DIR/$SECTOR/$ID
pwd
if [ -f "$SECTOR"_"$ID".asc ]
then
    echo "ERROR: found $(pwd)/${ID}.asc"
    echo "       AADG3 appears to have already finished successfully."
    exit 1
fi

AADG3 "$ID".in
