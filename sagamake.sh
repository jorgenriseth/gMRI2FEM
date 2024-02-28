#!/bin/bash
# Job name:
#SBATCH --job-name=snakerunner
#
# Project:
#SBATCH --account=nn9279k
#SBATCH --time='48:00:00'

# Max memory usage per task
#SBATCH --mem=500M

# Number of tasks (cores):
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 

# Logs
#SBATCH --output='logs/snakerunner_%j.out'

echo $@ 
echo "------"
echo "LOCALSCRATCH:" $LOCALSCRATCH
echo "SCRATCH": $SCRATCH 
snakemake $@ -p --profile "snakeprofiles/saga"  --configfile "snakeprofiles/saga/snakeconfig.yaml"