#! /bin/bash
## script to arrange all fastq files into batches for analysis. Combined files will be created for
## each barcode of interest for both passed and failed reads (e.g. pass_barcode01 or fail_barcode01).
## Samtools is used to quality-check and sort the file for each barcode, output is a sorted .bam-file.
## Sorted .bam-files are in addition converted to fasta-files. This serves as basis for ground truth-preparation 
## for each read and its respective base-modification.

# Get variables from command line
expName=$1 	# experiment name (date of exp usually)
bc=$2 		# barcode
refGenome=$3 	# full path to reference genome - passed from calling script as compund of $genomeFile/$referenceFiles
varSettingsFile=$4	# name of configuration file

# read Settings-file
source ${varSettingsFile}


####### modules to load ##########
#module load vital-it
#module load R/3.5.1
#module add UHTS/Analysis/minimap2/2.12;
#module add UHTS/Analysis/samtools/1.8;
############################################

#for bc in "${barcodesOfInterest[@]}" 
#do


#################################################
# Collecting reads from barcodes that were used #
#################################################
echo "collecting reads from folder of barcodes that were used..."

# merge all reads from particular barcode into single file (pass fail separately)
echo ${dataDir}/fast5files
#mkdir -p ${workDir}/qc/NanoStat

# need absolute paths for nanopolish index. get it from the summary file.
# summaryFile=${workDir}/fastqFiles/sequencing_summary.txt

if [ -d "${workDir}/barcoding/${bc}" ]; # check whether directory of specific barcode exists
then
    echo "passed reads..."
    cat ${workDir}/barcoding/${bc}/*.fastq > ${workDir}/barcoding/${expName}_pass_${bc}.fastq
    gzip ${workDir}/barcoding/${expName}_pass_${bc}.fastq

    #rm ${workDir}/bcFastq/${expName}_pass_${bc}.fastq
    #nanopolish index -s ${summaryFile} -d ${dataDir}/fast5Files ${workDir}/bcFastq/${expName}_pass_${bc}.fastq.gz
    #mkdir -p ${workDir}/qc/pass_${bc}
    #nanoQC ${workDir}/bcFastq/${expName}_pass_${bc}.fastq.gz -o ${workDir}/qc/pass_${bc}
    #NanoStat --fastq ${workDir}/bcFastq/${expName}_pass_${bc}.fastq.gz  --outdir ${workDir}/qc/NanoStat --name NanoStat_${expName}_pass_${bc}.txt --readtype 1D 
fi


############################
# aligning reads to genome #
############################

echo "aligning to genome using minimap2..."

mkdir -p ${workDir}/bamFiles

# map reads to genome with minimap2
minimap2 -ax map-ont $refGenome ${workDir}/barcoding/${expName}_pass_${bc}.fastq.gz | `# -a output in SAM-format; -ax map-ont for mapping of long noisy genomic reads, output in SAM-format` \
		samtools view -u -q 30 - | `# -u output uncompressed BAM (faster, as it will be piped anyway ("-" after -q option indicates redirection to standard output; -q MAPQ cutoff`\
		samtools sort -T pass_${bc} -o ${workDir}/bamFiles/${expName}_pass_${bc}.sorted.bam `sort SAM-file and output as BAM`

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
samtools fasta ${workDir}/bamFiles/${expName}_pass_${bc}.sorted.bam \
	> ${workDir}/bamFiles/${expName}_pass_${bc}.sorted.fasta
