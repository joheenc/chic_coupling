#!/bin/bash

#SBATCH --account=astro
#SBATCH -J leo_coupling

TLEFILE="starlink_2172.txt"
TSTART=59410 #15 Jul 2021
TSTOP=59590
LMAX=30

python3 pythonscripts/tle_to_orbit.py $TLEFILE $TSTART $TSTOP $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT
python3 pythonscripts/orbit_to_blm.py $SLURM_ARRAY_TASK_ID $LMAX
