#! /bin/bash

# Prepare mix of reads containing equal number of all specified barcodes
# Barcode selection (not necessarily all barcodes used in sequencing run)
# is specified in the varSettings.sh-file


varSettingsFile=$1 	# configuration file containing variable settings, paths etc.

############################################
# read settings and set up other variables #
############################################
source ${varSettingsFile}

bamDir=${workDir}/bamFiles
echo $bamDir

mkdir -p ${workDir}/barcodeMix_selectionForTraining
bcMixDir=${workDir}/barcodeMix_selectionForTraining
echo $bcMixDir

##################################################################################
# create summary file containing barcode, chromosome, leftmost mapping position, #
# length of read, MAPQ, and readID of ALL reads of ALL barcodes                  #
##################################################################################
# initiate summary file (tab-separated header only) #
echo -e "barcode_no"'\t'"chromosome"'\t'"leftmost_mapping_position"'\t'"readLength_bp"'\t'"MAPQ"'\t'"readID" \
	> ${bamDir}/readIDs_infoSummary_exp${expName}_allAlignedReads.txt

# continue summary file *
for i in ${barcodesStrippedAll[@]}
do
	samtools view ${bamDir}/${expName}_pass_barcode${i}.sorted.bam \
	| awk -v bcInd=$i -v OFS="\t" '{if (!($3=="chrM")) {print bcInd, $3, $4, length($10), $5, $1}}'
done >> ${bamDir}/readIDs_infoSummary_exp${expName}_allAlignedReads.txt
`# exclude reads from mitochondrial DNA`


###############################
# retrieve median read length #
###############################
readIDsummary=${bamDir}/readIDs_infoSummary_exp${expName}_allAlignedReads.txt
noReads=$(awk '{print $4}' ${readIDsummary} | tail -n +2 | wc -l)
median=$((noReads/2))
#median=$(( $(awk '{print $4}' ${readIDsummary} | tail -n +2 | wc -l) / 2 ))	# one-liner for median
medianLength=$(awk -v median=${median} '{print $4}' ${readIDsummary} | tail -n +2 | sort -g | tail -n ${median} | head -1)
echo -e "Experiment ${expName}\nnumber of aligned reads:\t${noReads}\nmedian read length [bp]:\t${medianLength}" | tee ${bamDir}/lengthInfo_AlignedReads.txt
# saving median length to file for use in subsequent scripts

##########################
# set read length to use #
##########################
if [ -z ${minReadLength} ]; then minReadLength=${medianLength}; fi


##############################################################################################
# retrieve minimum number available of reads with specific length per barcode and chromosome #
##############################################################################################
numReadsAvail=$( for i in ${barcodesStrippedSelection[@]}; \
		do \
		for j in ${chrOfInterest[@]}; \
			do \
			awk -v bc=$i -v chr=chr$j -v min=$minReadLength ' { if($4>=min && $1==bc && $2==chr) {print $0} }' ${bamDir}/readIDs_infoSummary_exp${expName}_allAlignedReads.txt \
			| cut -f 1,2 | sort | uniq -c ; \
			done; \
		done \
		| sort -g | head -1 | awk -v OFS='\t' ' {print $1, "(no. of reads)", "barcode" $2, $3}' | cut -f 1 )

echo "Minimum number of reads available per chromosome and barcode (minimum length ${minReadLength} bp): " ${numReadsAvail}


##################################
### set number of reads to use ###
##################################
if [ ${numReadsAvail} -lt ${minNumberOfReads} ]
then
    echo "Too few reads available"; echo "Requested: ${minNumberOfReads}, available ${numReadsAvail}"
    minNumberOfReads=$numReadsAvail  # set the minimal number of reads to the available number
    echo "Available number available used: ${minNumberOfReads}"
else
    echo "Sufficient reads available: ${minNumberOfReads} of ${numReadsAvail} used."
fi


#################################################################################################
# select the specific number of reads for each barcode and chromosome from individual BAM-files #
#################################################################################################
for i in ${barcodesStrippedSelection[@]}
do
	for j in ${chrOfInterest[@]}
	do
	samtools view -h ${bamDir}/${expName}_pass_barcode${i}.sorted.bam \
	| awk -v location=chr${j} -v min=${minReadLength} ' { if(NR<10) {print $0} else if (($3 == location) && (length($10) >= min)) {print $0} }' \
	| head -$((minNumberOfReads+9)) \
	| samtools view -b -h -o ${bamDir}/selection_exp${expName}_barcode${i}_chr${j}_min${minReadLength}bp.bam - 
	done
done
# comment on head -$((minNumberOfReads+9)): +9 serves to include the header lines


##############################
# merge individual BAM-files #
##############################
samtools merge ${bcMixDir}/selection_exp${expName}_barcodeMix_allChr_min${minReadLength}bp.bam  `# merged file` \
	${bamDir}/selection_exp${expName}_barcode*_chr*_min${minReadLength}bp.bam `# individual files to be merged` 


#############################################
# create summary file of the read selection #
#############################################
# initiate summary file (tab-separated header only) #
echo -e "barcode_no"'\t'"chromosome"'\t'"leftmost_mapping_position"'\t'"readLength_bp"'\t'"MAPQ"'\t'"readID" \
	> ${bcMixDir}/readIDs_infoSummary_exp${expName}_readSelection_allChr_min${minReadLength}bp.txt

for i in ${barcodesStrippedSelection[@]}
do
	for j in ${chrOfInterest[@]}
	do
		samtools view ${bamDir}/selection_exp${expName}_barcode${i}_chr${j}_min${minReadLength}bp.bam | \
		awk -v bcInd=$i -v location=chr${j} -v OFS="\t" '{ print bcInd, $3, $4, length($10), $5, $1 }'
	done >> ${bcMixDir}/readIDs_infoSummary_exp${expName}_readSelection_allChr_min${minReadLength}bp.txt
done


echo "Writing result files to " ${bamDir} " and " ${bcMixDir}

###########################################
# delete individual (temporary) BAM-files #
###########################################
rm ${bamDir}/selection_exp${expName}_barcode*_chr*_min${minReadLength}bp.bam  # delete individual files to save space
