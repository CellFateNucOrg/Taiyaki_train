#!/bin/bash

## Allocate resources
#SBATCH --time=0-04:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=32
#SBATCH --partition=all
##SBATCH --mail-user=
##SBATCH --mail-type=END,FAIL

## job name
#SBATCH --job-name="alignTBC"  # =alignTaiyakiBasecalls

# retrieve name of configuration file from command line
varSettingsFile=$1

# read Settings-file
source $varSettingsFile


# call script alignFastq for each barcode (=each entry $i of the job array)
./align_TaiyakiBasecalls.sh ${expName} $genomeFile/$referenceFiles $varSettingsFile
# location of reference to align to is a combination of folder (genomeFile) and actual file (referenceFiles)
	# if several modified reference files are used make use of this line:
	#./alignFastq_mod.sh ${expName} ${barcodesOfInterest[$i]} $genomeFile/$referenceFiles[$i]

