#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J read_depth


module load bioinfo-tools
module load samtools
module load BEDTools


# calculate read depth for coding sequence in genes


bam="$1"
regions=/home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_DEPTH/LsinapisSweM_chr_CDS_regions.txt
genes=/home/larshook/LarsH/FastZ/LsinapisSweM_gene_regions.bed

rm -f /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_DEPTH/"$1"-sorted.bam-unique.temp

while read line
do

  samtools depth \
	-r $line \
	-aa /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_MAPPING/MARK_DUPLICATES/"$1"-sorted.bam-unique.deduped.bam |\
		awk -v OFS="\t" '{print $1, $2-1, $2, $3}' >> /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_DEPTH/"$1"-sorted.bam-unique.temp
done <"$regions"


printf "chr pos depth start stop gene_id\n" > /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_DEPTH/"$1"-sorted.bam-unique.cov

bedtools intersect -a /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_DEPTH/"$1"-sorted.bam-unique.temp -b "$genes" -wa -wb |\
	awk '{print $1, $3, $4, $6, $7, $8}' |\
		sed 's/G/T/' >> /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_DEPTH/"$1"-sorted.bam-unique.cov


#rm -f /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_DEPTH/"$1"-sorted.bam-unique.temp
