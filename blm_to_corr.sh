#!/bin/bash

#SBATCH --account=astro
#SBATCH -J leo_coupling
#SBATCH --mem-per-cpu=120gb

LMAX=100
now=$(date)
echo "Start job at $now"

python3 pythonscripts/blm_to_coupling.py $LMAX
python3 pythonscripts/coupling_to_corr.py $LMAX

now=$(date)
echo "End job at $now"

