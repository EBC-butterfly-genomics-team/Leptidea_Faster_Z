#!/bin/bash -l

#SBATCH -A snic2022-5-34
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J BBmap_repair

module load bioinfo-tools
module load bbmap/38.61b
module load gzip

cd /home/larshook/LarsH/FastZ/Snps/Read_trimming
mkdir -p READS_FINAL

for sample in Swe-sin-101C_AGTTCC_L005 Swe-sin-102C_ATGTCA_L005 Swe-sin-1C_GCCAAT_L004 Swe-sin-2C_CAGATC_L004 Swe-sin-31C_ACTTGA_L004 Swe-sin-32C_GATCAG_L004 Swe-sin-61C_TAGCTT_L004 Swe-sin-62C_GGCTAC_L005 Swe-sin-91C_CTTGTA_L005 Swe-sin-92C_AGTCAA_L005
do

  cp /home/larshook/LarsH/FastZ/Snps/Read_trimming/FILTERED_READS/"$sample"_R*_001_val_*.tagged_filter.fastq.gz $SNIC_TMP

  cd $SNIC_TMP

  repair.sh -Xmx6G \
        in1="$sample"_R1_001_val_1.tagged_filter.fastq.gz \
	in2="$sample"_R2_001_val_2.tagged_filter.fastq.gz \
	out1="$sample"_R1_final.fq \
	out2="$sample"_R2_final.fq \
	repair


  gzip *final.fq

  cp *final.fq.gz /home/larshook/LarsH/FastZ/Snps/Read_trimming/READS_FINAL

done
