#!/bin/bash
## Script to basecall from fast5-files using the model generated with Taiyaki
## Model: Last available checkpoint from traing round 2
## Output: Basecalls and information about modified bases

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
# activate taiyaki #
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

mkdir -p ${workDir}/taiyakiBasecalls_CpG_GpC_modelTest

# Taiyaki basecalling
basecall.py --device cuda --modified_base_output ${workDir}/taiyakiBasecalls_CpG_GpC_modelTest/basecalls_modeltest_CpG_GpC_${expName}.hdf5 `# --device cuda --> run on GPUs;  output of modified base information` \
	${dataDir}/fast5Files `# directory containing raw data (fast5-files` \
	${workDir}/training2_run1001_CpG_GpC_modified_bases/model_checkpoint_00017.checkpoint `# trained model file (final checkpoint from 2nd round of training`\
	> ${workDir}/taiyakiBasecalls_CpG_GpC_modelTest/basecalls_modeltest_CpG_GpC_${expName}.fa `# output basecalls redirected (by > ) to this file`

