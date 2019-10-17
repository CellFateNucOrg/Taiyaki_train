#! /bin/bash
## script to sort fastq-Files by barcode
## .bashrc should include location of bin directory inside Guppy directory and export it as ${GUPPYDIR}
## output will be in barcoding-subfolder of working directory


# retrieve name of configuration file from command line
varSettingsFile=$1

# read Settings-file
source ${varSettingsFile}



######################
# sorting by barcode #
######################

guppy_barcoder -i ${workDir}/fastqFiles -s ${workDir}/barcoding/ --barcode_kits EXP-NBD104 --recursive
