#!/usr/bin/env bash

# when using this command, you must export the SECTOR variable
# and specify the array range
# e.g. sbatch --export=SECTOR="north" --array=0-99 MESA_sector_array.slurm

#SBATCH --nodes 1
#SBATCH --tasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 4G
#SBATCH --time 00:40:00
#SBATCH --output slurm/MESA_sector_array_%A_%a.out
#SBATCH --error slurm/MESA_sector_array_%A_%a.err

module purge; module load bluebear  # this line is required
module load scipy/0.19.0-python-3.5.2

export OMP_NUM_THREADS=2
echo OMP_NUM_THREADS $OMP_NUM_THREADS
echo $SECTOR
export MESA_DIR=/rds/projects/2017/ballwh-stellar-modelling/mesa/7385
cd $SLURM_SUBMIT_DIR/$SECTOR/$(printf "%05d" $SLURM_ARRAY_TASK_ID)
pwd
if [ -f final.profile.GYRE ]
then
    echo "ERROR: found $(pwd)/final.profile.GYRE"
    echo "       MESA appears to have already finished successfully."
    exit 1
fi

if ! [ -L star ]
then
    ln -s /rds/projects/2017/ballwh-stellar-modelling/s4tess/mesa_work/star
fi

# if MESA log exists, then see which try we were on, and retry
# otherwise, this is the first try

MAXTRIES=95
if [ -f mesa.log ]
then
    TRIES=$(head -n 1 mesa.log | awk '{print $2}')
else
    TRIES=1
fi

echo TRY $TRIES > mesa.log
./rn >> mesa.log
# while grep --quiet convergence mesa.log && ((TRIES <= MAXTRIES))
while [ ! -f final.profile.GYRE ] && ((TRIES <= MAXTRIES))
do
    python3 ../../../../scripts/modify_Tc.py inlist_run 1000 -d
    ((TRIES++))
    echo TRY $TRIES > mesa.log
    ./rn >> mesa.log
done
