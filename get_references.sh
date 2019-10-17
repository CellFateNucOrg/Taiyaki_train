#!/bin/bash

# Task: generate references (Taiyaki-module get_refs_from_sam.py)
# Taiyaki-module extracts a snippet of the reference for each mapped read
# which serve as read-specific references

###################################
# Get variables from command line #
###################################
varSettingsFile=$1 	# configuration file containing variable settings, paths etc.


###################################
# read and set up other variables #
###################################
source ${varSettingsFile}
# e.g. $workDir, $dataDir

# set working directories for this step
bamDir=${workDir}/bamFiles
bcMixDir=${workDir}/barcodeMix_selectionForTraining
groundTruthDir=${workDir}/groundTruth_barcodeMix_allChr

if [ -z ${minReadLength} ]; then minReadLength=$(tail -1 ${bamDir}/lengthInfo_AlignedReads.txt | cut -f 2); fi  # use median read length (determined in barcodeMix_preparation.sh) if a minimum read length is not specified in varSettings.sh
# use median read length (determined in barcodeMix_preparation.sh) if a minimum read length is not specified in varSettings.sh
echo "Starting extraction of read-specific references..."
echo "minimum read length used [bp]: " ${minReadLength}

# create directory for ground truth
mkdir -p ${workDir}/groundTruth_barcodeMix_allChr


outFile=${groundTruthDir}/selection_exp${expName}_barcodeMix_allChr_min${minReadLength}bp_taiyaki_referenceMatches.fasta  # read-specific references


####################
# activate taiyaki #
####################
source ${TAIYAKI_DIR}/venv/bin/activate


####################################
# generate read-specific reference #
####################################
echo "Retrieving BAM-file containing read selection from" ${bcMixDir}
### call Taiyaki module (get_refs_from_sam.py) to generate read-specific references ###
${TAIYAKI_DIR}/bin/get_refs_from_sam.py \
	$genomeFile/$referenceFiles  `# genome reference file` \
	${bcMixDir}/selection_exp${expName}_barcodeMix_allChr_min${minReadLength}bp.bam  `# BAM-file containing the reads selected for training (generated by barcodeMix_preparation.sh). The Taiyaki module get_refs_from_sam.py chooses read-specific references for each of the reads in this file.` \
	--output ${outFile} `# file containing a snippet of the reference for each mapped read, i.e. a read-specific reference ` 

# Comment: The get_refs_from_sam.py module will find read-specific reference for most (but NOT all) reads
echo "Saving results (Fasta-file) to" ${outFile}

