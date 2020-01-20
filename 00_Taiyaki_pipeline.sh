#!/bin/bash

#  There are several options for the '--dependency' flag that depend on the status of Job1. e.g.
# --dependency=afterany:Job1	Job2 will start after Job1 completes with any exit status
# --dependency=after:Job1	Job2 will start any time after Job1 starts
# --dependency=afterok:Job1	Job2 will run only if Job1 completed with an exit status of 0
# --dependency=afternotok:Job1	Job2 will run only if Job1 completed with a non-zero exit status

# Example usage: sbatch 00_pipeline_final.sh -v varSettings.sh -e 4 -x 6  to start from step 4 (barcodeMix preparation) to step 6 (ground truth preparation)
# or sbatch 00_pipeline_final.sh -v varSettings.sh -p -c -e 1 -x 7 for start from basecalling (entryPoint -e 1) and proceed up to merging of ground truth files (exitPoint -x 7), including creation of scaling parameters (-c) and pycoQC quality check (-p)

print_usage() {
  printf "Usage:	-v [path/name of variable settings file]\n"		# name of configuration file from command line
  printf "	-c flag triggering creation of scaling parameters (default no)\n"	# set flag to enable creation of scaling parameters
  printf "	-p flag triggering pycoQC (default no)\n"		# set flag to launch pycoQC
  printf "	-e [number from 1 through 8] entry point to pipeline\n"
  printf "	-x [number from 1 through 8] exit point from pipeline\n"
  printf "		steps:	1 basecalling \n"
  printf "			2 barcode separation (demultiplexing) \n"
  printf "			3 alignment \n"
  printf "			4 barcode mix generation for training \n"
  printf "			5 reference sequence retrieval \n"
  printf "			6 generation of ground truth \n"
  printf "			7 merging ground truth of individual barcodes \n"
  printf "			8 pre-training (read mapping) \n"
#  printf "	-m flag triggering first model training \n"
#  printf "	-n flag triggering second model training \n"
}

if [ $# -eq 0 ]; then
    echo "No arguments provided"
    print_usage; exit 1
fi

####################
### Set defaults ###
####################
createScalingParams=false
runPycoQC=false
training_1=false
training_2=false

while getopts 'v:cpe:x:mn' flag; do
  case "${flag}" in
	v) varSettingsFile="${OPTARG}" ;;
    c) createScalingParams=true ;;
    p) runPycoQC=true ;;
    e) entryPoint="${OPTARG}" ;;
    x) exitPoint="${OPTARG}" ;;
    m) training_1=true ;;
    n) training_2=true ;;
    *) print_usage
       exit 1 ;;
  esac
done


######################
# read settings-file #
######################
source ${varSettingsFile}

###########################################################
### Ready for pycoQC? (successful basecalling required) ###
###########################################################

if [ -f "${workDir}/fastqFiles/sequencing_summary.txt" ];  # checks whether pycoQC can be started right away
then pycoQC_readyToGo="true"; else pycoQC_readyToGo="false";
fi


###########################################################################
### Processing steps of the Taiyaki pipeline (up to training the model) ###
###########################################################################
steps=( 01_runBasecall.sh 01b_runBarCode.sh 02_runAlignFastq.sh 03_barcodeMix.sh 03a_getReferences.sh \
		03b_groundTruthReference.sh 03c_mergeGroundTruthFiles.sh 05_taiyaki_prepare_mapped_reads.sh )
# 04_create_scaling_params.sh and pycoQC are not included in the array as they are done in parallel


#################
# echo settings #
#################
echo "Taiyaki workflow pipeline started..."
echo -e "variable settings file:" '\t' $varSettingsFile  # -e enables the interpretation of backslash escapes
echo -e "working directory:" '\t' $workDir
echo -e "data directory:" '\t' $dataDir
echo -e "create scaling parameters:" $createScalingParams '\t' "run pycoQC:" $runPycoQC
echo -e "entry point:" $entryPoint "(${steps[$((entryPoint-1))]})" '\t' "exit point:" $exitPoint "(${steps[$((exitPoint-1))]})"


############################################################################################################
# start with generation of per read-scaling parameters in parallel (takes long and requires raw data only) #
############################################################################################################
if [[ $createScalingParams == true ]] && [ -z "$(ls -A ${workDir}/per_read_scaling_parameters 2> /dev/null)" ]  # -z checks for zero length; the second operator hence checks whether directory exists and/or is empty (true if it does not exist; true if it exists but is empty; false if it exists and is not empty
then
	echo "Generating per read-scaling parameters..."
	echo "Writing files to ${workDir}/per_read_scaling_parameters..."
	sbatch 04_create_scaling_params.sh ${varSettingsFile}
else
	echo "Directory per_read_scaling_parameters/ already exists and is not empty (or parameter creation disabled)."
	echo "Check directory contents before re-running 04_create_scaling_params.sh separately!"
fi


#############################################################
# declare array variable jobIDs used to store SLURM job IDs #
#############################################################
declare -a jobIDs


###################################################
# Taiyaki pipeline - jobs depend on previous ones #
###################################################

if !( [ -z "$entryPoint" ] && [ -z "$exitPoint" ] )  # entry and exit points need to be defined
then
	if [[ "$entryPoint" == "$exitPoint" ]] # just a single step is carried out, there is no dependencies
	then
		echo "Single step submitted: ${steps[$((entryPoint-1))]}"
		sbatch ${steps[$((entryPoint-1))]} ${varSettingsFile}
	elif [[ $entryPoint < $exitPoint ]]  # more than one step (plus, order of entry and exit points is correct)
	then
		for (( i = ${entryPoint}; i <= ${exitPoint}; i++))  # loop over all steps from entryPoint to exitPoint
		do
			if [ "$i" == "$entryPoint" ]  # the entry step - starts independently
			then
				eval "jobIDs[$((i-1))]"="$(sbatch ${steps[$((i-1))]} ${varSettingsFile} | cut -d ' ' -f 4)"
				sleep 1
				echo "First of ${#steps[@]} steps: ${steps[$((i-1))]}  [ SLURM jobID ${jobIDs[$((i-1))]} ]"
							# The SLURM jobID is assigned to an array variable in the array jobIDs
							# The cut command is required to extract the jobID from the return message ("Submitted batch job [jobID]")
			else  # steps following the entry step - start depends on successful completion of previous job
				eval "jobIDs[$((i-1))]"="$(sbatch --dependency=afterok:${jobIDs[-1]} ${steps[$((i-1))]} ${varSettingsFile} | cut -d ' ' -f 4)"
				sleep 1
				echo "Step $(($i-$entryPoint+1)) of ${entryPoint} to ${exitPoint} submitted (from a total of ${#steps[@]} steps): ${steps[$((i-1))]} - [ SLURM jobID ${jobIDs[$((i-1))]} ]"
							# The next job is started only after successful completion of the previous one (which has the SLURM ID stored on the last position, i.e. index [-1], of the jobIDs array)
			fi
		done
	fi
else
	echo "No entry / exit points defined."
fi


##############################################################################################
# quality control (pycoQC) - starts once sequencing_summary is available (after basecalling) #
##############################################################################################
if ([ "$runPycoQC" == "true" ] && [ "$pycoQC_readyToGo" == "true" ])
then
	sbatch 01a_pycoQC.sh ${varSettingsFile}
	echo "Starting pycoQC.."
	echo "Saving quality report to ${workDir}/fastqFiles..."
elif [ "$runPycoQC" == "true" ] && [ "$pycoQC_readyToGo" == "false" ]  # the pycoQC_readyToGo-flag is false
then
	echo "Waiting for job ${jobIDs[$((entryPoint-1))]} to finish..."
	sbatch --dependency=afterok:${jobIDs[$((entryPoint-1))]} 01a_pycoQC.sh ${varSettingsFile}  # basecalling is the step required for pycoQC, so if the first step is completed pycoQC can start
	echo "Starting pycoQC.."
	echo "Saving quality report to ${workDir}/fastqFiles..."
fi


##################
# Model training #
##################
# check if training is requested #
if [ "${training_1}" == "true" ] || [ "${training_2}" == "true" ]; then doTraining=true; else exit; fi


if !( [ -z "$entryPoint" ] && [ -z "$exitPoint" ] )
then
	let lastPipelineJob=${#JobIDs[*]}-1
	waitingRequired=true
#	dependency='--dependency=afterok:${jobIDs[${lastPipelineJob}]}'
else
	waitingRequired=false
fi

#echo doTraining ${doTraining}
#echo lastPipelineJob ${lastPipelineJob}
#echo waitingRequired ${waitingRequired}
#echo training_1 ${training_1}
#echo training_2 ${training_2}


declare -a jobID_train
declare -a varSettingsArray
varSettingsArray=( ${varSettingsFile} )
echo varArray ${varSettingsArray[0]}


if [ "${doTraining}" == "true" ] && [ ${waitingRequired} == "true" ]
then
	if [ "${training_1}" == "true" ] && [ "${training_2}" =0 "true" ]
	then
		eval "jobID_train[1]"="$(sbatch --job-name=train_1 --dependency=afterok:${jobIDs[${lastPipelineJob}]} taiyaki_modelTraining.sh  ${varSettingsFile}  1)"
#		echo "$(sbatch --job-name=train_1 ${dependency} taiyaki_modelTraining.sh  ${varSettingsFile}  1)"
		echo "${jobID_train[1]}"
		sleep 1
		eval "jobID_train[2]"="$(sbatch --job-name=train_2 --dependency=afterany:${jobID_train[1]}]}  taiyaki_modelTraining.sh  ${varSettingsFile}  2)"
		#echo "jobID_train_2"="$(sbatch --job-name=train_2 --dependency=afterany:${jobID_train_1}]}  taiyaki_modelTraining.sh  ${varSettingsFile}  2)"
	elif [ "${training_1}" == "true"  ]
	then
		eval "jobID_train[1]"="$(sbatch --job-name=train_1 --dependency=afterok:${jobIDs[${lastPipelineJob}]} taiyaki_modelTraining.sh  ${varSettingsFile}  1)"
		sleep 1
	else
		eval "jobID_train_2"="$(sbatch --job-name=train_2  taiyaki_modelTraining.sh  ${varSettingsFile}  2)"
	fi
elif [ "${doTraining}" == "true" ] && [ ${waitingRequired} == "false" ]
then
	if [ "${training_1}" == "true" ] && [ "${training_2}" = "true" ]
	then
		eval "jobID_train[1]"="$(sbatch --job-name=train_1  taiyaki_modelTraining.sh  ${varSettingsFile}  1)"
#		echo "$(sbatch --job-name=train_1  taiyaki_modelTraining.sh  ${varSettingsFile}  1)"
		echo "train1 JobID: ${jobID_train[1]}"
		sleep 2
		eval "jobID_train[2]"="$(sbatch --job-name=train_2 --dependency=afterany:${jobID_train[1]}]}  taiyaki_modelTraining.sh  ${varSettingsFile}  2)"
#		echo "jobID_train[2]"="$(sbatch --job-name=train_2 --dependency=afterany:${jobID_train_1}]}  taiyaki_modelTraining.sh  ${varSettingsFile}  2)"
	elif [ "${training_1}" == "true"  ]
	then
		eval "jobID_train_1"="$(sbatch --job-name=train_1  taiyaki_modelTraining.sh  ${varSettingsFile}  1)"
		sleep 1
	else
		eval "jobID_train_2"="$(sbatch --job-name=train_2  taiyaki_modelTraining.sh  "${varSettingsArray[0]}"  2)"
	fi
fi

sleep 3

