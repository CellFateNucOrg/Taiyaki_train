#! /bin/bash
## script to arrange all fastq files into batches for analysis. Combined files will be created for
## each barcode of interest for both passed and failed reads (e.g. pass_barcode01 or fail_barcode01).
## Samtools is used to quality-check and sort the file for each barcode, output is a sorted .bam-file.
## Sorted .bam-files are in addition converted to fasta-files. This serves as basis for ground truth-preparation 
## for each read and its respective base-modification.

# Get variables from command line
expName=$1 	# experiment name (date of exp usually)
refGenome=$2 	# full path to reference genome - passed from calling script as compund of $genomeFile/$referenceFiles
varSettingsFile=$3	# name of configuration file

# read Settings-file
source ${varSettingsFile}

############################
# update $workDir-variable #
############################
workDir=${workDir%/*}/modelData_and_training_barcodeMix/taiyakiBasecalls_CpG_GpC_modelTest
echo "Taiyaki basecalls in ${workDir}"



####### modules to load ##########
#module load vital-it
#module load R/3.5.1
#module add UHTS/Analysis/minimap2/2.12;
#module add UHTS/Analysis/samtools/1.8;
############################################


############################
# aligning reads to genome #
############################

echo "Aligning Taiyaki basecalls to genome using minimap2..."

mkdir -p ${workDir}/bamFiles_alignment_Taiyaki_basecalls_to_reference

# map reads to genome with minimap2
minimap2 -ax map-ont $refGenome  ${workDir}/basecalls_modeltest_CpG_GpC_${expName}.fa | `# -a output in SAM-format; -ax map-ont for mapping of long noisy genomic reads, output in SAM-format` \
		samtools sort -o ${workDir}/bamFiles_alignment_Taiyaki_basecalls_to_reference/${expName}_Taiyaki_basecalls_CpG_GpC_alignment.sorted.bam `sort SAM-file and output as BAM`

# samtools view -q: MAPQ is the MAPping Quality. It equals −10 log10 Pr{mapping position is wrong}, rounded to the nearest integer. A value 255 indicates that the mapping quality is not available.
# If the probability of a correct match is 0.999, the MAPQ score is 30.


# filter reads with flag=2308: unmapped (4) + secondary alignment (256) + supplementary alignment (2048) [the latter category is the main problem]
#minimap2 -ax map-ont $refGenome ${workDir}/bcFastq/${expName}_fail_${bc}.fastq.gz | samtools view -F 2308 -b | samtools sort -T fail_${bc} -o ${workDir}/bamFiles/${expName}_fail_${bc}.sorted.bam

#echo "index bam file ..."
#samtools index ${workDir}/bamFiles/${expName}_pass_${bc}.sorted.bam
#samtools index ${workDir}/bamFiles/${expName}_fail_${bc}.sorted.bam



##################################################
# creating fasta files for all barcode-bam-files #
##################################################

# convert sorted .bam files to fasta-format (for read-specific reference creation in case of base modifications)
samtools fasta ${workDir}/bamFiles_alignment_Taiyaki_basecalls_to_reference/${expName}_Taiyaki_basecalls_CpG_GpC_alignment.sorted.bam \
	> ${workDir}/bamFiles_alignment_Taiyaki_basecalls_to_reference/${expName}_Taiyaki_basecalls_CpG_GpC_alignment.sorted.fasta
