#! /bin/bash
## Taiyaki model training

## Allocate resources
#SBATCH --time=7-00:00:00
#SBATCH --partition=all
#SBATCH --gres=gpu:1 ## reserve GPU (1 means yes, i.e. GPU will be reserved)
#SBATCH --mem=96GB

## job name
##SBATCH --job-name="train1stModel"


# retrieve variables from command line
varSettingsFile=$1 	# name of configuration file
trainingRound=$2


# run model training
./taiyaki_modelTraining.sh ${varSettingsFile} ${trainingRound}
