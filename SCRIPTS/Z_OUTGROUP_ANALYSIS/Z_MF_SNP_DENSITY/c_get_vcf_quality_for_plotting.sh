#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J vcf_quality


ml bioinfo-tools
ml bcftools


cd /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/SNPS/MPILEUP

for sample in P14502_103 P14502_104 P14458_107_S28 P14458_108_S29 P14458_103_S24 P14458_104_S25
do
  echo "QUAL DP" > "$sample"_quality.txt
  bcftools view $sample.bcf |\
	grep -v "#" |\
		sed 's/DP=//;s/;/ /g' |\
			awk '{print $6, $8}' >> "$sample"_quality.txt
done
