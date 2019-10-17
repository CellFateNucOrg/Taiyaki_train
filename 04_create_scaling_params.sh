#!/bin/bash
## script to create per-read scaling parameters using Taiyaki module generate_per_read_params.py

## Allocate resources
#SBATCH --time=7-0:00:00
#SBATCH --partition=all
#SBATCH --exclude=izblisbon,izbcotonou
#SBATCH --mem=96G
#SBATCH --cpus-per-task=32

## job name
#SBATCH --job-name="createParams"

# retrieve configuration file from command line
varSettingsFile=$1 	# name of configuration file

# read in the run-specific settings
source ${varSettingsFile}


####################
# activate taiyaki #
####################
source ${TAIYAKI_DIR}/venv/bin/activate


mkdir -p ${workDir}/per_read_scaling_parameters

######################################
# Create Per-read Scaling Parameters #
######################################

python ${TAIYAKI_DIR}/bin/generate_per_read_params.py --jobs 32 ${dataDir}/fast5Files \
	> ${workDir}/per_read_scaling_parameters/read_params_${expName}.tsv
# --jobs [integer] sets the number of threads used, here 32 - this should match with --cpus-per-task set in the allocation 
