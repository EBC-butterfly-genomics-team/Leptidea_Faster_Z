#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J vcf_filter

module load bioinfo-tools
module load bcftools
module load python/3.8.7


cd /home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/SNPS/MPILEUP


for sample in P14502_103 P14502_104 P14458_107_S28 P14458_108_S29 P14458_103_S24 P14458_104_S25 P14502_105 P14502_106
do
  bcftools filter \
        -i 'GT="het" && QUAL>30 && sum(DP4)>25 && sum(DP4)<200 && (DP4[0]+DP4[1])>10 && (DP4[2]+DP4[3])>10' \
        -o "$sample"-filtered_200.vcf.gz \
        -Oz \
	"$sample".bcf
done
