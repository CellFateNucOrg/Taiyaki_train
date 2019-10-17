#!/bin/bash

## Allocate resources

## Allocate resources
#SBATCH --time=0-2:00:00

################################
# Difference to other scripts! #
################################

#SBATCH --array=1-4			## Submit n jobs (#SBATCH --array=1-n), with n being the number of barcodes in barcodesSelection,
#####################################################################################################################
# Do NOT use the number of barcodes in barcodesOfInterest as before but the number of barcodes in barcodesSelection #
#####################################################################################################################

#SBATCH --partition=all
#SBATCH --mem=16G
#SBATCH --cpus-per-task=24


## job name
#SBATCH --job-name="groundTruth"

# retrieve configuration file from command line
varSettingsFile=$1 	# name of configuration file

# read in the run-specific settings
source $varSettingsFile


let i=$SLURM_ARRAY_TASK_ID-1	# TASK_ID-1 because indices (e.g. of barcode array) start at 0
# =index in barcode array, used to determine type of modification in readModifications.sh


# call script readModifications.sh for each barcode (=each entry $i of the job array)
./readModifications.sh  ${i}  $varSettingsFile

# Depending on the barcode used, specific bases in the sequence will be modified, generating a ground truth-reference
# Which barcode translates into which modification has to be aligned between varSettings.sh and ground_truth_and_mods.sh

