#! /bin/bash
## script to basecall all fast5 in the folder called fast5files and its recursive directories with Guppy
## .bashrc should include location of bin directory inside Guppy directory and export it as ${GUPPYDIR}
## uses barcoding
## output will be in ${workDir}/fastqFiles, the fastqFiles subfolder of the working directory

# retrieve configuration file from command line
varSettingsFile=$1 	# configuration file containing variable settings, paths etc.

# load configuration
source ${varSettingsFile}


#######################
# basecalling (guppy) #
#######################

guppy_basecaller --input_path ${dataDir}/fast5Files --save_path ${workDir}/fastqFiles \
		--records_per_fastq 200000 --recursive  --qscore_filtering --min_qscore 3  \
		--gpu_runners_per_device 5 --cpu_threads_per_caller 2 \
		--device cuda:0 --config ${guppyConfigFile}
