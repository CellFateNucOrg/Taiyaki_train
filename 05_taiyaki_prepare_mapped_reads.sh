#! /bin/bash
## script to mapped reads (alignment between raw signal and reference sequence)

## Allocate resources
#SBATCH --time=3-00:00:00  # time format: days-hours:minutes:seconds
#SBATCH --partition=gpu
#SBATCH --mem=96G
#SBATCH --cpus-per-task=32

## job name
#SBATCH --job-name="mapReads"

#################################################
# retrieve configuration file from command line #
#################################################
varSettingsFile=$1 	# name of configuration file

######################################
# read settings and set up variables #
######################################
source ${varSettingsFile}
if [ -z ${minReadLength} ]; then minReadLength=$(tail -1 ${bamDir}/lengthInfo_AlignedReads.txt | cut -f 2); fi  # use median read length (determined in barcodeMix_preparation.sh) if a minimum read length is not specified in varSettings.sh

groundTruthDir=${workDir}/groundTruth_barcodeMix_allChr
jointGroundTruthFile=selection_exp${expName}_barcodeMix_allChr_min${minReadLength}bp_taiyaki_groundTruth.fasta

mkdir -p ${workDir}/modelData_and_training
modelDir=${workDir}/modelData_and_training

####################
# activate taiyaki #
####################
source ${TAIYAKI_DIR}/venv/bin/activate


###########################
# Create Mapped Read File #
###########################
echo 'Mapping reads (alignment between raw signal and reference sequence...'
echo 'Creating mapped reads for mixed barcode-input...'
echo 'Running on CPUs...'

modNamesSelection=$(for i in ${indicesOfSelection[*]}; \
	do \
		printf "%s__" ${modificationsOfInterest[$i]}; \
	done)	# prints specific elements of array modificationsOfInterest (specified in array indicesOfSelection) \
			# and uses __ as delimiter for concatenation

modNamesSelection=$(echo ${modNamesSelection%__})  # %_ strips the trailing "__" from the variable name

# specify alternative base identifiers using --mod altBaseName base baseMod, e.g. --mod Z C 5mC  --mod Y A 6mA
# (one --mod-expression per modified base)

${TAIYAKI_DIR}/bin/prepare_mapped_reads.py --device cpu --recursive --no-overwrite --jobs 32 --mod Z C 5mC  \
	${dataDir}/fast5Files `# location of fast5 files of this particular run (or possibly a subset thereof)` \
	${workDir}/per_read_scaling_parameters/read_params_${expName}.tsv  `# file containing the read-specific scaling parameters, from generate_per_read_params.py`\
	${modelDir}/modbases_${expName}_ModificationsUsedForTraining_${modNamesSelection}.hdf5  `# output file (mapped reads)`\
	${TAIYAKI_DIR}/models/mGru_flipflop_remapping_model_r9_DNA.checkpoint  `# pre-trained model / entry model` \
	${groundTruthDir}/selection_exp${expName}_barcodeMix_allChr_min${minReadLength}bp_taiyaki_groundTruth.fasta `# fasta file containing a reference specific for each read marked up with modified base information` \


