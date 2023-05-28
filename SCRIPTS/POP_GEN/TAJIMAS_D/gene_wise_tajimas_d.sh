#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J tajimas_d


# calulate gene wise Tajima's D using gene coordinates and gene specific CDS snps...


ml bioinfo-tools
ml vcftools
ml BEDTools


cd /home/larshook/LarsH/FastZ/TAJIMAS_D
mkdir -p OUT
cd OUT


bedfile=/home/larshook/LarsH/FastZ/PI/LsinapisSweM.all.maker.genes.bed
vcf=/home/larshook/LarsH/FastZ/TAJIMAS_D/CDS_all_sites_LsinapisSweM-Swe-sin_final_all-filtered_snps.vcf.gz

# loop through bedfile with gene coords
# take snps in each gene with bedtools and make temp vcf
# calculate Tajima's D for snps in each gene

while read chr start end gene
do 
  grep -w "#CHROM" $vcf > $gene.vcf
  echo $chr $start $end $gene |\
	awk -v OFS="\t" '{print $1, $2, $3, $4}' |\
		bedtools intersect \
			-a /home/larshook/LarsH/FastZ/TAJIMAS_D/CDS_all_sites_LsinapisSweM-Swe-sin_final_all-filtered_snps.vcf.gz \
			-b - >> $gene.vcf
  vcftools \
  	--vcf ${gene}.vcf \
  	--out ${gene} \
  	--TajimaD ${end} \
	--chr ${chr}

done <$bedfile


rm -f *.vcf
rm -f *.log
rm -f ../gene_wise_Tajimas_D.txt


# concatenate results into one file...

for file in LSSWEG*.Tajima.D
do 
  gene=$(echo $file | sed 's/.Tajima.D//')
  TajD=$(grep -v CHROM $file | awk '{print $4}')
  echo $gene $TajD >> ../gene_wise_Tajimas_D.txt
done

rm -f LSSWEG*.Tajima.D
