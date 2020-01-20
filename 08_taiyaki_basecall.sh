#!/bin/bash
## Script to basecall from fast5-files using the model generated with Taiyaki
## Model: Last available checkpoint from traing round 2
## Output: Basecalls (fasta-file) and information about modified bases (hdf5-file)

## Allocate resources
#SBATCH --time=2-0:00:00
#SBATCH --partition=all
#SBATCH --gres=gpu:1 ## reserve GPU (1 means yes, i.e. GPU will be reserved)
#SBATCH --mem=128GB

## job name
#SBATCH --job-name="basecTaiyaki"

# get name of Settings-file from command line
varSettingsFile=$1	# name of configuration file

# read Settings-file
source ${varSettingsFile}

####################
# activate Taiyaki #
####################
source ${TAIYAKI_DIR}/venv/bin/activate


#######################
# Taiyaki basecalling #
#######################
# Taiyaki had performed two rounds of training: 
# the first round down-weights learning the modified bases in favour a good canonical call, 
# the second round then focuses on learning the conditional prediction of whether a base is modified.

# Now, using the model generated in the second training, basecalling is performed on the raw data (fast5 files).
# The basecalls produced use the canonical base alphabet, information about putative modifed base calls
# is written out to the specified .hdf5-file


# name of training directory and model files#
modNamesSelection=$(for i in ${indicesOfSelection[*]}; \
	do printf "%s__" ${modificationsOfInterest[$i]}; \
	done)	# concatenates specific elements of array modificationsOfInterest (specified in array indicesOfSelection) with "__" as delimiter for concatenation
modNamesSelection=$(echo ${modNamesSelection%__})  # %__ strips the trailing "__" from the variable name

### directory containing model file (checkpoint files of model training) ###
training_2_Dir=${modelDir}/training2_exp${expName}_${modNamesSelection}_${modelDirNameAppendix[1]}
### Checkpoint file to be used (result of training round 2) ###
finalCheckpoint=$(ls ${training_2_Dir}/*.checkpoint | tail -1)

### input directory (fast5-files)
# directory containing raw data (fast5-files) for basecalling ###
# path of inputDirTaiyakiFast5 is sourced from varSettings.sh-file (hence commented out below in line 54), 
# but could be set below as well. Change directory accordingly if different set of files is to be called
# inputDirTaiyakiFast5=${dataDir}/fast5Files

### output directory
taiyakiBasecallDir=${workDir}/Taiyaki_basecalls_exp${expName}_${modNamesSelection}

mkdir -p ${taiyakiBasecallDir}

#######################
# Taiyaki basecalling #
#######################
basecall.py --device cuda --modified_base_output ${taiyakiBasecallDir}/Taiyaki_modifiedBases_exp${expName}_${modNamesSelection}.hdf5 `# --device cuda --> run on GPUs;  output of modified base information (hdf5-file)` \
	${inputDirTaiyakiFast5} `# directory containing raw data (fast5-files), not necessarily the files used for training` \
	${finalCheckpoint} `# trained model file (usually final checkpoint from 2nd round of training)`\
	> ${taiyakiBasecallDir}/Taiyaki_canonicalBasecalls_${expName}_${modNamesSelection}.fa `# output basecalls (fasta-file) redirected (by > ) to this file`


