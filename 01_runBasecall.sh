#!/bin/bash

## Allocate resources
#SBATCH --time=10:00:00
#SBATCH --partition=all
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32		# would only apply if run on CPUs

## job name
#SBATCH --job-name="dSMFguppy"

# retrieve configuration file from command line
varSettingsFile=$1 	# name of configuration file


./basecallGuppy.sh $varSettingsFile
