#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J make_alternate_cds


# pick out cds sequence from assembly fasta based on gff coordinates for reference, variant and snp masked assembly...

module load bioinfo-tools
module load BEDTools


gff=/proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.CDS.gff
main_path=/home/larshook/LarsH/FastZ
reference_path=/proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies

ref_fasta=LsinapisSweM
snp_fasta=LsinapisSweM_variant
masked_fasta=LsinapisSweM_snp_masked

# reference...

sed 's/ID=/ /;s/-RA/ /' $gff |\
       awk -v OFS="\t" '{print $1, $2, $9, $4, $5, $6, $7, $8}' |\
		bedtools getfasta \
			-fi $reference_path/"$ref_fasta".fasta \
			-bed - \
			-s \
			-name |\
				sed s'/:/ /' |\
					awk '{print $1}' |\
						fold -w 60 > $main_path/"$ref_fasta"_CDS.fasta

awk '$1 ~ /^>/ {if (a[$1]++ < 1) print $1}; $1 !~ /^>/ {print}' $main_path/"$ref_fasta"_CDS.fasta > $main_path/"$ref_fasta"_concatenated_CDS.fasta


# variant...

sed 's/ID=/ /;s/-RA/ /' $gff |\
	awk -v OFS="\t" '{print $1, $2, $9, $4, $5, $6, $7, $8}' |\
		bedtools getfasta \
       		-fi $main_path/"$snp_fasta".fasta \
			-bed - \
       			-s \
			-name |\
               			sed s'/:/ /' |\
	                      		awk '{print $1}' |\
       	                      			fold -w 60 > $main_path/"$snp_fasta"_CDS.fasta

awk '$1 ~ /^>/ {if (a[$1]++ < 1) print $1}; $1 !~ /^>/ {print}' $main_path/"$snp_fasta"_CDS.fasta > $main_path/"$snp_fasta"_concatenated_CDS.fasta


# snp masked...

sed 's/ID=/ /;s/-RA/ /' $gff |\
	awk -v OFS="\t" '{print $1, $2, $9, $4, $5, $6, $7, $8}' |\
		bedtools getfasta \
	        	-fi $main_path/"$masked_fasta".fasta \
			-bed - \
        		-s \
			-name |\
                		sed s'/:/ /' |\
                        		awk '{print $1}' |\
                                		fold -w 60 > $main_path/"$masked_fasta"_CDS.fasta

awk '$1 ~ /^>/ {if (a[$1]++ < 1) print $1}; $1 !~ /^>/ {print}' $main_path/"$masked_fasta"_CDS.fasta > $main_path/"$masked_fasta"_concatenated_CDS.fasta
