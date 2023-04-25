#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH -J read_depth


module load bioinfo-tools
module load samtools

bam="$1"
regions=LsinapisSweF_chr_CDS_regions.txt

printf "chr pos depth cds_midpoint\n" > /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_DEPTH/"$1"-sorted.bam-unique.cov

while read line
do

  samtools depth \
	-r $line \
	-aa /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_MAPPING/MARK_DUPLICATES/"$1"-sorted.bam-unique.deduped.bam |\
	awk -v id=$line '{print $1, $2, $3, id}' |\
		sed 's/:/ /;s/-/ /' |\
			awk -v OFMT='%f' '{print $1, $2, $3, ($5+$6)/2}' >> /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_DEPTH/"$1"-sorted.bam-unique.cov
done </home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_DEPTH/"$regions"
