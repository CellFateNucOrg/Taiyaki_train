#!/bin/bash

# Task: generate references (Taiyaki-module get_refs_from_sam.py)
# Taiyaki-module extracts a snippet of the reference for each mapped read

## Allocate resources
#SBATCH --time=0-5:00:00
#SBATCH --partition=gpu  # redirects to izbdelhi without automatically blocking GPU resources
##SBATCH --exclude izblisbon,izbcotonou # currently, Taiyaki is only installed on izbdelhi, so excluding izblisbon and izbcotonou will guide the script there
## Do not use SBATCH --gres=gpu:1 to be directed to izbdelhi. No GPUs are required in this step. Use --exclude option (previous line) instead
#SBATCH --mem=96G
#SBATCH --cpus-per-task=32

## job name
#SBATCH --job-name="getRefs"

# retrieve configuration file from command line
varSettingsFile=$1 	# name of configuration file

# read the settings
source $varSettingsFile

####################################################################################################
# call script to generate read-specific references (snippet of the reference for each mapped read) #
####################################################################################################
./get_references.sh $varSettingsFile

