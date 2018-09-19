#!/usr/bin/env bash

# when using this command, you must export the RUN variable
# and specify the array range
# e.g. sbatch --export=RUN="alpha" --array=0-99 GYRE_sector_array.slurm

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 1G
#SBATCH --time 00:45:00
#SBATCH --output slurm/GYRE_sector_array_%A_%a.out
#SBATCH --error slurm/GYRE_sector_array_%A_%a.err
#SBATCH --account=ballwh-stellar-modelling

module purge; module load bluebear  # this line is required

ID=$(($START + $SLURM_ARRAY_TASK_ID))
echo "ID               $ID"

export OMP_NUM_THREADS=1
echo "OMP_NUM_THREADS  $OMP_NUM_THREADS"
echo "RUN              $RUN"

TMPFILE=$(mktemp /scratch/${USER}_MESA_${SLURM_JOBID}_${ID}.XXXXXX)

export GYRE_DIR=/rds/homes/b/ballwh/code/gyre/fork
cd $SLURM_SUBMIT_DIR/$RUN/$(printf "%05d" $ID)
pwd
if [ -f gyre_summary.txt ]
then
    echo "ERROR: found $(pwd)/gyre_summary.txt"
    echo "       GYRE appears to have already finished successfully."
    exit 1
fi

$GYRE_DIR/bin/gyre gyre.in > $TMPFILE
mv -v $TMPFILE gyre.log
