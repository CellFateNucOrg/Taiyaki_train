#!/bin/bash

# show/plot training progress of Taiyaki training

## Allocate resources
#SBATCH --time=00:10:00
#SBATCH --partition=all
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

## job name
#SBATCH --job-name="plotTrain"

# retrieve configuration file from command line
varSettingsFile=$1 	# name of configuration file
trainingRound=$2


#################
# read settings #
#################
source $varSettingsFile


# directory to retrieve training results from
modNamesSelection=$(for i in ${indicesOfSelection[*]}; \
	do printf "%s__" ${modificationsOfInterest[$i]}; \
	done)	# concatenates specific elements of array modificationsOfInterest (specified in array indicesOfSelection) with "__" as delimiter for concatenation
modNamesSelection=$(echo ${modNamesSelection%__})  # %__ strips the trailing "__" from the variable name

trainingIndex=$((trainingRound-1))
trainingDir=${modelDir}/training${trainingRound}_exp${expName}_${modNamesSelection}_${modelDirNameAppendix[${trainingIndex}]}

# name of progress file (.png-image)
timeStamp=$(date +%Y-%m-%d_%H%M%S)
progressFile=training_${trainingRound}_progress_${timeStamp}


#usage: plot_training.py [-h] [--mav MAV] [--upper_y_limit UPPER_Y_LIMIT]
#                        [--lower_y_limit LOWER_Y_LIMIT]
#                        output folder plus file name  input_directories [input_directories ...]


####################
# activate Taiyaki #
####################
source ${TAIYAKI_DIR}/venv/bin/activate


# create progress plot
${TAIYAKI_DIR}/misc/plot_training.py --lower_y_limit 0.0001 --mav 1  `# lower y-limit set to 0.0001 (almost 0); moving average (mav) to smooth line, choose value for mav between 1 and approx. 10, or leave out --mav option`\
			${trainingDir}/${progressFile} `# path and name of progress file` \
			${trainingDir} `# input directory (containing model.all_loss.txt and model.log` 


