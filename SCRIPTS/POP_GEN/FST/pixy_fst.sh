#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 100:00:00
#SBATCH -J pixy


module load bioinfo-tools
module load pixy
module load samtools


bedfile=LsinapisSweM_w100kb_s10kb.bed
fold="$1"

cd /home/larshook/LarsH/FastZ
mkdir -p FST
cd FST
mkdir -p BED


cp /home/larshook/LarsH/FastZ/PI/LsinapisSweM_"$fold"-fold_sites.txt $SNIC_TMP  
cp /home/larshook/SCRIPTS/FastZ/FST/population_file.txt $SNIC_TMP
cp /home/larshook/LarsH/FastZ/"$bedfile" $SNIC_TMP
cp /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2/LsinapisSweM-sin-juv_final_all-filtered_q30.vcf.gz $SNIC_TMP

cd $SNIC_TMP

tabix LsinapisSweM-sin-juv_final_all-filtered_q30.vcf.gz

pixy \
      	--stats fst dxy pi \
        --vcf LsinapisSweM-sin-juv_final_all-filtered_q30.vcf.gz \
        --populations population_file.txt \
	--sites_file LsinapisSweM_"$fold"-fold_sites.txt \
	--bed_file "$bedfile" \
        --n_cores 2

cp -r * /home/larshook/LarsH/FastZ/FST/"$fold"_FOLD/BED
