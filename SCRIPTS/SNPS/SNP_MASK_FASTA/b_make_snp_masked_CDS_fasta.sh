#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J make_alternate_cds


# use snp masked and variant genome fasta to generate snp masked and variant CDS sequences using annotation gff coordiantes

module load bioinfo-tools
module load BEDTools


gff=/proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.CDS.gff
main_path=/home/larshook/LarsH/FastZ



snp_fasta=LsinapisSweM_variant
masked_fasta=LsinapisSweM_snp_masked

sed 's/ID=/ /;s/-RA/ /' $gff |\
	awk -v OFS="\t" '{print $1, $2, $9, $4, $5, $6, $7, $8}' |\
		bedtools getfasta \
			-fi $main_path/$snp_fasta.fasta \
			-bed - \
			-s \
			-name |\
				sed s'/:/ /' |\
					awk '{print $1}' |\
						fold -w 60 > $main_path/"$snp_fasta"_CDS.fasta


sed 's/ID=/ /;s/-RA/ /' $gff |\
        awk -v OFS="\t" '{print $1, $2, $9, $4, $5, $6, $7, $8}' |\
		bedtools getfasta \
        		-fi $main_path/$masked_fasta.fasta \
        		-bed - \
        		-s \
			-name |\
				sed s'/:/ /' |\
					awk '{print $1}' |\
						fold -w 60 > $main_path/"$masked_fasta"_CDS.fasta
