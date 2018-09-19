#!/usr/bin/env bash

# when using this command, you must export the RUN variable
# and specify the array range
# e.g. sbatch --export=RUN="alpha" --array=0-99 AADG3_sector_array.slurm

#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mem 256M
#SBATCH --time 00:10:00
#SBATCH --output slurm/postprocess_sector_array_%A_%a.out
#SBATCH --error slurm/postprocess_sector_array_%A_%a.err
#SBATCH --account=ballwh-stellar-modelling

module purge; module load bluebear  # this line is required
module load scipy/0.19.0-python-3.5.2

# export /rds/homes/b/ballwh/python/
ID=$(($START + $SLURM_ARRAY_TASK_ID))

echo "ID              $ID"
echo "RUN             $RUN"
ID5=$(printf "%05d" $ID)
cd $SLURM_SUBMIT_DIR
python3 ../scripts/add_white_noise.py $RUN/$ID5
python3 ../scripts/sectorize.py $RUN/$ID5
ls $RUN/$ID5/"$ID5"_WN_*.asc | xargs -t -n1 python3 ../scripts/asc_to_fits.py
