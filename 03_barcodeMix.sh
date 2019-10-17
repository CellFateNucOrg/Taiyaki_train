#!/bin/bash

## Allocate resources

## Allocate resources
#SBATCH --time=0-1:00:00
#SBATCH --partition=all
##SBATCH --exclude izblisbon,izbcotonou # currently, Taiyaki is only installed on izbdelhi, so excluding izblisbon and izbcotonou will guide the script there
## Do not use SBATCH --gres=gpu:1 to be directed to izbdelhi. No GPUs are required in this step. Use --exclude option (previous line) instead
#SBATCH --mem=32G
##SBATCH --cpus-per-task=4
## job name
#SBATCH --job-name="bcMix"

# retrieve configuration file from command line
varSettingsFile=$1 	# name of configuration file

##########################
# read variable settings #
##########################
source $varSettingsFile


#################################################
# call script for preparation of the barcodeMix #
#################################################
./barcodeMix_preparation.sh $varSettingsFile

