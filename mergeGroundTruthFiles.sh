#!/bin/bash

# Merge the individual ground truth fasta-files of all barcodes used in the training selection
# see indicesOfSelection in varSettings-file for the barcodes used for the training selection


###################################
# Get variables from command line #
###################################
varSettingsFile=$1 	# configuration file containing variable settings, paths etc.


###################################
# read and set up other variables #
###################################
source ${varSettingsFile}
if [ -z ${minReadLength} ]; then minReadLength=$(tail -1 ${bamDir}/lengthInfo_AlignedReads.txt | cut -f 2); fi  # use median read length (determined in barcodeMix_preparation.sh) if a minimum read length is not specified in varSettings.sh


groundTruthDir=${workDir}/groundTruth_barcodeMix_allChr
indivGroundTruthFiles=selection_exp${expName}_allChr_min${minReadLength}bp_taiyaki_groundTruth_*.fasta  # ground truth: read-specific references using alternative letters for the modified bases
jointGroundTruthFile=selection_exp${expName}_barcodeMix_allChr_min${minReadLength}bp_taiyaki_groundTruth.fasta


#############################################
# concatenate individual ground truth files #
#############################################
echo "Generating single ground truth-file..."
cat ${groundTruthDir}/${indivGroundTruthFiles} > ${groundTruthDir}/${jointGroundTruthFile} &&
echo "Done. Ground truth file for Taiyaki model training saved in ${groundTruthDir}."
