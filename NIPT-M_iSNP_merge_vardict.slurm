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
# Update: 2023.5.8
# Usage: Merge GVCF files and call variants
# Note: Apply for gatk and vardict output vcfs
# *******************************************


# ******************************************
# 0. Setup
# ******************************************
# load all the GVCF files
GVCF_1="/dssg/home/acct-medkwf/medkwf4/results/NIPT-M/Lib20230429-01/output.g.vcf.gz"
GVCF_2="/dssg/home/acct-medkwf/medkwf4/results/NIPT-M/Lib20230429-02/output.g.vcf.gz"
GVCF_3="/dssg/home/acct-medkwf/medkwf4/results/NIPT-M/Lib20230429-03/Lib20230429-03.vcf"

# General configs
NT=16 #NT=$(nproc) #number of threads to use in computation, set to number of cores in the server
PWD="/dssg/home/acct-medkwf/medkwf4/results/NIPT-M"
START_DIR="$PWD/merge_vardict" #Determine where the output files will be stored

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

# ******************************************
# 4. Merge GCVF
# ******************************************
picard MergeVcfs \
I=${GVCF_1} \
I=${GVCF_2} \
I=${GVCF_3} \
O=combined.raw.vcf
