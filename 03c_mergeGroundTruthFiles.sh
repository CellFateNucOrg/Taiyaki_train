#!/bin/bash

## Allocate resources

## Allocate resources
#SBATCH --time=0-2:00:00
#SBATCH --partition=all
#SBATCH --mem=96G
#SBATCH --cpus-per-task=32


## job name
#SBATCH --job-name="mergeGT"

# retrieve configuration file from command line
varSettingsFile=$1 	# name of configuration file

# read in the run-specific settings
source $varSettingsFile

# call script readModifications.sh for each barcode (=each entry $i of the job array)
./mergeGroundTruthFiles.sh  $varSettingsFile

