#!/bin/bash

#######################
# guppy configuration #
#######################
guppyConfigFile=/home/jsemple/mysoftware/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg


#################################
# directory of Taiyaki software #
#################################
# may be available as global variable already #
#TAIYAKI_DIR=/home/jsemple/mysoftware/taiyaki


######################
# name of experiment #
######################
expName="20190830"


################################
# data and working directories #
################################
dataDir=/home/jsemple/myimaging/Jenny/pipeline_test
workDir=/home/jsemple/myimaging/Jenny/pipeline_test

bamDir=${workDir}/bamFiles
bcMixDir=${workDir}/barcodeMix_selectionForTraining
groundTruthDir=${workDir}/groundTruth_barcodeMix_allChr
modelDir=${workDir}/modelData_and_training


#####################################
# location of reference genome file #
#####################################
# location used is a combination of folder (genomeFile) and actual file (referenceFiles),
# i.e. $genomeFile/$referenceFiles
# used in scripts 02_runAlignFastq.sh and 05_taiyaki_prepare_mapped_reads.sh
genomeFile=/mnt/imaging.data/pmeister/ce11
referenceFiles=genome_ce11_no_linebreaks.fa


######################################################
# barcodes and modifications used in this experiment #
#####################################################################################################################
# Fill in carfully:																									#
#		chrOfInterest			--> chromosomes the training set shall be compiled from								#
#		barcodesOfInterest		--> barcodes used in experiment														#
#		modificationsOfInterest --> names of the modifications corresponding to each barcode						#
#		numberOfModifications	--> number of modifications for each barcode group									#
#		indicesOfSelection		--> indices of the entries in array barcodesOfInterest to be used for training set	#
#		motifs and corresponding modifications for each barcode group												#
#		minReadLength			--> minimum length of reads to be used for training set (defaults to median length)	#
#		minNumberOfReads		--> minimum number of reads to be used per chromosome and barcode					#
#									(if not reached lowest number of available reads per chromosome and barcode		#
#									is used																			#
#####################################################################################################################

### chromosomes of interest ###
chrOfInterest=( I II III IV V X )	# chromosomes of interest (specified e.g. to ignore reads from mitochondrial DNA or focus on a specific chromosome)

### barcodes used ###
barcodesOfInterest=( barcode01 barcode02 barcode03 barcode04 barcode06 )
barcodesStrippedAll=( $(echo ${barcodesOfInterest[@]#barcode}) )  # strips leading string "barcode" from the barcodes to be used for training
		#barcodesPlain=( $(echo ${barcodesStripped[@]#0}) )  # strips leading string "0" if applicable (sequential stripping ensures that "barcode" is stripped in any case) - barcodesPlain was not used in the end as iteration works with the stripped barcodes as well

# identifiers (names) of the modifications corresponding to the barcodes #
modificationsOfInterest=( "no_mod" "CpG" "GpC" "CpG_GpC" "6mA" )  # Array containing the names of the modification motifs corresponding to the barcodes

# the number of modifications in each barcode group, i.e. the number of methyltransferases used
numberOfModifications=( 0 1 1 2 1 )

### motifs and modifications corresponding to the barcodes ###
### motif01 / 02 are the motifs recognized by the first and (potentially) second methyltransferase used during library preparation, respectivels ###
# Likewise, modification01 / 02 correspond to the modified motifs (motif01 / 02) using alternative letters, respectively, e.g. using Z instead of methylated C
# The position in array motif or modification must match the respective position in array barcodesOfInterest
		motif01=( ' '  'CG'  'GC'  'CG'  'A' )	# array of motifs to be modified
 modification01=( ' '  'ZG'  'GZ'  'ZG'  'Y' )	# array of identifiers for the modified motifs using alternative letters for the modified base

		motif02=( ' '  ' '   ' '   'GC'  ' ' )	# array of motifs to be modified
 modification02=( ' '  ' '   ' '   'GZ'  ' ' )	# array of identifiers for the modified motifs using alternative letters for the modified base
# The second array may be mainly empty as only in few cases a second methyltransferase is used.
# The entries containing just a space ' ' are nevertheless required for the substitution to work properly!


##################################
# specifications of training set #
##################################
### Training selection ###
# Indices corresponding to the respective entries in the array barcodesOfInterest
# These indices are used to generate the "selection"-variables
# Be aware that indix numbering start at 0, so e.g. index 3 stands for the element in position number 4) or 0, 2, 3 to use elements 1, 3 and 4
indicesOfSelection=( 0 1 2 3 )

minReadLength=500  # minimum length of reads selected for training set
	# if minReadLength is not specified (or line is commented out), the median read length will be used
minNumberOfReads=3000  # desired number of reads per chromosome and barcode to be combined as training set
	# if number available is exceeded by minNumberOfReads, the number available will be used (so putting an arbitrarilly high number will lead to use of all available reads)

### barcodesSelection ###
for i in ${indicesOfSelection[*]}; do barcodesSelection+=("${barcodesOfInterest[$i]}"); done
for i in ${indicesOfSelection[*]}; do barcodesStrippedSelection+=("${barcodesStrippedAll[$i]}"); done

# identifiers (names) of the selection #
for i in ${indicesOfSelection[*]}; do modificationsOfInterestSelection+=("${modificationsOfInterest[$i]}"); done


#############################
# model training parameters #
#############################
trainingModFactor=( 0.01 1.0 )  # --mod_factor argument controls the proportion of the training loss attributed to the modified base output stream in comparison to the canonical base output stream
modelDirNameAppendix=( "_canonicalCall" "_modifiedBaseRecognition" )


########## don't edit below this line ###############
# note: the .bashrc should point to the location of nanopolish with the variable
# NANOPOLISH_DIR

export NUMBC=${#barcodesOfInterest[@]}
