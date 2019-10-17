#!/bin/bash
## script for pycoQC quality control after basecalling
## output will be in ${workDir}/fastqFiles, i.e. subfolder fastqFiles of the working directory

# retrieve configuration file from command line
varSettingsFile=$1 	# configuration file containing variable settings, paths etc.

# load configuration
source ${varSettingsFile}



####################
#  Quality control #
####################

pycoQC -f ${workDir}/fastqFiles/sequencing_summary.txt -o ${workDir}/fastqFiles/pycoQC_${expName}_output.html
