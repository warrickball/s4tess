#!/usr/bin/env bash

# when using this command, you must export the RUN variable
# and specify the array range
# e.g. sbatch --export=RUN="alpha" --array=0-99 AADG3_sector_array.slurm

#SBATCH --nodes 1
#SBATCH --tasks 1
#SBATCH --mem 1G
#SBATCH --time 01:00:00
#SBATCH --output slurm/AADG3_sector_array_%A_%a.out
#SBATCH --error slurm/AADG3_sector_array_%A_%a.err

module purge; module load bluebear  # this line is required

echo $RUN
ID=$(printf "%05d" $SLURM_ARRAY_TASK_ID)
cd $SLURM_SUBMIT_DIR/$RUN/$ID
pwd
if [ -f "$ID".asc ]
then
    echo "ERROR: found $(pwd)/${ID}.asc"
    echo "       AADG3 appears to have already finished successfully."
    exit 1
fi

AADG3 "$ID".in
