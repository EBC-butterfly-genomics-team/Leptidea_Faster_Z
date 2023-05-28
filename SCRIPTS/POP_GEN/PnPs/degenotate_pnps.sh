#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J P-sites

ml bioinfo-tools
ml samtools


# use degenotate to get Pn and Ps...
# output file mk.tsv contains estimates


# remove "*" segmental deletion annotation as these makes degenotate crash... 

zcat /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2/LsinapisSweM-sin-juv_final_all-snps-filtered_q30.vcf.gz |\
	 sed 's/,\*//;s/\*,//' > /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2/LsinapisSweM-sin-juv_final_all-snps-filtered_q30_del_trimmed.vcf

bgzip /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2/LsinapisSweM-sin-juv_final_all-snps-filtered_q30_del_trimmed.vcf

tabix /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2/LsinapisSweM-sin-juv_final_all-snps-filtered_q30_del_trimmed.vcf.gz

/home/larshook/LarsH/SOFTWARE/Degenotate/Degenotate/bin/degenotate.py
	-a /home/larshook/LarsH/FastZ/DN_DS/LsinapisSweM.all.maker.genes_rebuilt.gff \
	-g /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/LsinapisSweM.fasta \
	-v /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2/LsinapisSweM-sin-juv_final_all-snps-filtered_q30_del_trimmed.vcf.gz \
	-u /home/larshook/SCRIPTS/FastZ/PN_PS/outgroup.txt \
	-o /home/larshook/LarsH/FastZ/PN_PS \
	--overwrite
