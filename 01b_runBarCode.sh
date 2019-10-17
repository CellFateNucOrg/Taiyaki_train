#!/bin/bash

## Allocate resources
#SBATCH --time=7-00:00:00
#SBATCH --partition=all
#SBATCH --exclude izblisbon,izbcotonou  # to run on izbdelhi without occupying GPUs
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32

## job name
#SBATCH --job-name="barcodeGuppy"

# retrieve name of configuration file from command line
varSettingsFile=$1

./barcodeGuppy.sh $varSettingsFile
