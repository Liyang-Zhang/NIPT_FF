#!/bin/bash

#SBATCH --mail-type=end
#SBATCH --mail-user=zly18810602991@163.com
#SBATCH --job-name=test
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -J NIPT-M_iSNP
#SBATCH --ntasks-per-node=8
#SBATCH --output=%j.o
#SBATCH --error=%j.e

# *******************************************
# Author: Liyang.Zhang
# Version: v1.0.0
# Update: 2023.5.6
# Usage: Script to perform DNA seq variant calling
#        using an target region sample with fastq files
#        named 1.fastq.gz and 2.fastq.gz
# *******************************************

# ******************************************
# 0. Setup
# ******************************************
# For each submission:
# Update with the fullpath location of your sample fastq
# Change the following setting for ctDNA-EQA:
#     1. SM
#     2. FASTQ_FOLDER
#     3. FASTQ_1
#     4. FASTQ_2

# Update with the fullpath location of your sample fastq
SM="Lib20230429-03" #sample name
RGID="rg_$SM" #read group ID
PL="ILLUMINA" #or other sequencing platform
FASTQ_FOLDER="/dssg/home/acct-medkwf/medkwf9/data/cgdata/DNA/reprod/NIPT-M/YUEER-20230429-L-01-2023-05-021545/Sample_JZ23020675-20230429-Lib20230429-03"
FASTQ_1="$FASTQ_FOLDER/JZ23020675-20230429-Lib20230429-03_combined_R1.fastq.gz"
FASTQ_2="$FASTQ_FOLDER/JZ23020675-20230429-Lib20230429-03_combined_R2.fastq.gz" #If using Illumina paired data

# General configs
NT=16 #NT=$(nproc) #number of threads to use in computation, set to number of cores in the server
PWD="/dssg/home/acct-medkwf/medkwf4/results/NIPT-M"
START_DIR="$PWD/$SM" #Determine where the output files will be stored

# Change to the working directory
WORKDIR="$START_DIR"
mkdir -p $WORKDIR
cd $WORKDIR

# actiavte conda environment
module load miniconda3
source activate dna

# Update with the location of the reference data files
FASTA_DIR="/dssg/home/acct-medkwf/medkwf4/database/GATK/hg19"
FASTA="$FASTA_DIR/ucsc.hg19.fasta"
KNOWN_DBSNP="$FASTA_DIR/dbsnp_138.hg19.vcf.gz"
KNOWN_INDELS="$FASTA_DIR/1000G_phase1.indels.hg19.sites.vcf.gz"
KNOWN_MILLS="$FASTA_DIR/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
INTERVAL_FILE="/dssg/home/acct-medkwf/medkwf4/script/NIPT-M_iSNP/region_file/NIPT-M_iSNP.bed"

#Vardict directory
vardict_path="/dssg/home/acct-medkwf/medkwf4/.conda/envs/R_test/share/vardict-2018.09.21-0"

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000

#bwa mem -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL" -K 100000000 -t $NT $FASTA $FASTQ_1 $FASTQ_2 | samtools view -1 - -o align.bam
#samtools sort align.bam -o sort.bam
#samtools index sort.bam

# ******************************************
# 2. Mark duplicates
# ******************************************
#picard MarkDuplicates \
#INPUT=sort.bam \
#OUTPUT=dedup.bam \
#METRICS_FILE=dedup.metrics \
#VALIDATION_STRINGENCY=SILENT \
#OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
#ASSUME_SORT_ORDER="queryname" \
#CREATE_MD5_FILE=false
#
#samtools index dedup.bam

# ******************************************
# 3. Calling SNPs (HaplotypeCaller)
# ******************************************
#gatk --java-options "-Xmx128G" HaplotypeCaller \
#-R ${FASTA} \
#-I dedup.bam \
#-O output.g.vcf.gz \
#-L ${INTERVAL_FILE} \
#-ERC GVCF \

# ******************************************
# 4. Calling SNPs (Vardict)
# ******************************************
conda deactivate
source activate R_test
vardict -G ${FASTA} \
-f 0.001 \
-N ${SM} \
-b dedup.bam \
-c 1 \
-S 2 \
-E 3 \
-g 4 \
${INTERVAL_FILE} | ${vardict_path}/teststrandbias.R | ${vardict_path}/var2vcf_valid.pl -N ${SM} -E -f 0.001 > ${SM}.vcf

