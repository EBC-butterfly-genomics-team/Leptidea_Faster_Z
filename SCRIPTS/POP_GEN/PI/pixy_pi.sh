#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 50:00:00
#SBATCH -J pixy

module load bioinfo-tools
module load pixy
module load samtools
module load BEDTools


fold=N
window=100000
bedfile=LsinapisSweM_w100kb_s10kb.bed

cp /home/larshook/LarsH/FastZ/"$bedfile" $SNIC_TMP
cp /home/larshook/SCRIPTS/FastZ/PI/population_file.txt $SNIC_TMP
cp /home/larshook/LarsH/FastZ/PI/LsinapisSweM_"$fold"-fold_sites.txt $SNIC_TMP
cp /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2/LsinapisSweM-Swe-sin_final_all-filtered_q30.vcf.gz $SNIC_TMP

cd $SNIC_TMP


tabix LsinapisSweM-Swe-sin_final_all-filtered_q30.vcf.gz

pixy \
      	--stats pi \
        --vcf LsinapisSweM-Swe-sin_final_all-filtered_q30.vcf.gz \
        --populations population_file.txt \
        --sites_file LsinapisSweM_"$fold"-fold_sites.txt \
        --bed_file "$bedfile" \
        --n_cores 1

mv pixy_pi.txt pixy_pi_"$fold"-fold_"$bedfile"_q30.txt

cp pixy_pi_"$fold"-fold_"$bedfile"_q30.txt /home/larshook/LarsH/FastZ/PI/
