#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH -J vcf_quality


ml bioinfo-tools
ml bcftools


cd /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller/round2
sample=LsinapisSweM-sin-juv_final_all-filtered_q30


#cd /home/larshook/LarsH/FastZ/Snps/HaplotypeCaller
#sample=LsinapisSweM-Ire-juv_final_snps



# print all different info fields...

echo "quality score" > "$sample"_snp_quality.txt

for i in AF DP FS MQ QD
do
  bcftools view -v snps "$sample".vcf.gz |\
  	grep -v "##" |\
		awk 'NR>1 {print $8}' |\
			sed 's/;/ /g' |\
	 			awk -v OFS="\n" '{$1=$1}1' |\
	 				sed 's/=/ /;s/,/ /' |\
	 					awk '{if (NF==4) print $1, $2"\n"$1, $3"\n"$1, $4; else if (NF==3) print $1, $2"\n"$1, $3; else if (NF==2) print $1, $2}' |\
	 						grep -w "$i" >> "$sample"_snp_quality.txt
done
