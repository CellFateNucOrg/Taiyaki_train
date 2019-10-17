#!/bin/bash

# Taiyaki had been used to select read-specific references (get_refs_from_sam.py)
# These reference sequences are modified now by introducing barcode-specific base modifications
# to generate the ground truth for the model training


###################################
# Get variables from command line #
###################################
selectionIndex=$1	# index position in barcode array (already adjusted to start at 0 in the calling script)
varSettingsFile=$2 	# configuration file containing variable settings, paths etc.


###################################
# read and set up other variables #
###################################
source ${varSettingsFile}
if [ -z ${minReadLength} ]; then minReadLength=$(tail -1 ${bamDir}/lengthInfo_AlignedReads.txt | cut -f 2); fi  # use median read length (determined in barcodeMix_preparation.sh) if a minimum read length is not specified in varSettings.sh


######################################################
# barcodes and modifications used in this experiment #
######################################################
modName=${modificationsOfInterest[$selectionIndex]}
modNumber=${numberOfModifications[$selectionIndex]}
bcName=${barcodesStrippedAll[$selectionIndex]}

bcMixDir=${workDir}/barcodeMix_selectionForTraining
groundTruthDir=${workDir}/groundTruth_barcodeMix_allChr

readIDfile=${bcMixDir}/readIDs_infoSummary_exp${expName}_readSelection_allChr_min${minReadLength}bp.txt
taiyakiRefFile=${groundTruthDir}/selection_exp${expName}_barcodeMix_allChr_min${minReadLength}bp_taiyaki_referenceMatches.fasta  # read-specific references, generated by Taiyaki module get_refs_from_sam.py
groundTruthFile=${groundTruthDir}/selection_exp${expName}_allChr_min${minReadLength}bp_taiyaki_groundTruth_bc${bcName}_${modName}.fasta  # ground truth: read-specific references using alternative letters for the modified bases



######################
# base modifications #
######################
echo "Replacing modified bases with alternative letters in fasta-file..."
echo "barcode: " ${barcodesStrippedAll[$selectionIndex]}

printf "%-3s %-8s %-4s\n" ${motif01[$selectionIndex]} "exchanged with " ${modification01[$selectionIndex]}
if [ ${modNumber} == 2 ]; then
	printf "%-3s %-8s %-4s\n" ${motif02[$selectionIndex]} "exchanged with " ${modification02[$selectionIndex]}
fi

echo -n "command used: " sed s/${motif01[$selectionIndex]}/${modification01[$selectionIndex]}/g
if [ ${modNumber} == 2 ]; then
	echo " |" sed s/${motif02[$selectionIndex]}/${modification02[$selectionIndex]}/g
fi
echo

### introducing barcode-specific modifications ###
for ID in $(grep ^${barcodesStrippedAll[$selectionIndex]} ${readIDfile} | cut -f 6)
do
	grep -A1 ${ID} ${taiyakiRefFile} \
	| sed "s/${motif01[$selectionIndex]}/${modification01[$selectionIndex]}/g" \
	| sed "s/${motif02[$selectionIndex]}/${modification02[$selectionIndex]}/g"
done > ${groundTruthFile}

