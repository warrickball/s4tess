#!/usr/bin/env bash

# when using this command, you must export the RUN variable
# and specify the array range
# e.g. sbatch --export=RUN="alpha" --array=0-99 AADG3_sector_array.slurm

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 1G
#SBATCH --time 02:00:00
#SBATCH --output slurm/AADG3_sector_array_%A_%a.out
#SBATCH --error slurm/AADG3_sector_array_%A_%a.err
#SBATCH --account=ballwh-stellar-modelling

module purge; module load bluebear  # this line is required

ID=$(($START + $SLURM_ARRAY_TASK_ID))
echo "ID               $ID"

echo "RUN              $RUN"

TMPFILE=$(mktemp /scratch/${USER}_MESA_${SLURM_JOBID}_${ID}.XXXXXX)

ID5=$(printf "%05d" $ID)
cd $SLURM_SUBMIT_DIR/$RUN/$ID5
pwd
if [ -f "$ID5".asc ]
then
    echo "ERROR: found $(pwd)/${ID5}.asc"
    echo "       AADG3 appears to have already finished successfully."
    exit 1
fi

AADG3 "$ID5".in -v > $TMPFILE
mv -v $TMPFILE AADG3.log
