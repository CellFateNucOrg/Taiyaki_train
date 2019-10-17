#!/bin/bash

## Allocate resources
#SBATCH --time=3-00:00:00
#SBATCH --partition=all
#SBATCH --mem=256G
#SBATCH --cpus-per-task=32

## job name
#SBATCH --job-name="pycoQC"

# retrieve name of configuration file from command line
varSettingsFile=$1

./pycoQC.sh $varSettingsFile
