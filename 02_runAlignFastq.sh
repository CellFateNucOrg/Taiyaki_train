#!/bin/bash

## Allocate resources
#SBATCH --time=0-04:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=32
#SBATCH --array=1-5
#SBATCH --partition=all
##SBATCH --mail-user=
##SBATCH --mail-type=END,FAIL
## you should submit as many jobs as there are barcodes in barcodesOfInterest
## (don't forget to include unclassfied in barcodesOfInterst in the varSettings.sh file) - not done in previous scripts?

## job name
#SBATCH --job-name="npAlign"

# retrieve name of configuration file from command line
varSettingsFile=$1

# read Settings-file
source $varSettingsFile


let i=$SLURM_ARRAY_TASK_ID-1
#-1 because the indices are starting at 0

# align fastq-files using minimap2
# if several reference files are used: array of barcodesOfInterest and reference files needs to be aligned
# 	(e.g. barcode matching CmG must match the CmG-modified reference), position indicated
# 	by variable $i (only if an array of modified references is used)
# location of reference to align to is a combination of folder (genomeFile) and actual file (referenceFiles)

# call script alignFastq for each barcode (=each entry $i of the job array)
./alignFastq.sh ${expName} ${barcodesOfInterest[$i]} $genomeFile/$referenceFiles $varSettingsFile
	# if several modified reference files are used make use of this line:
	#./alignFastq_mod.sh ${expName} ${barcodesOfInterest[$i]} $genomeFile/$referenceFiles[$i]


# if barcode is one with a spikein, align also to other genomes
#for bc in ${bcWithSpikeIn[@]}
#do
#    if [ "${barcodesOfInterest[$i]}" == "$bc" ]; then
#        ./alignFastq_spikeIn.sh ${expName} ${barcodesOfInterest[$i]} ${lambdaFile}
#        ./alignFastq_spikeIn.sh ${expName} ${barcodesOfInterest[$i]} ${phiXfile}
#    	exit
#    fi
#done
