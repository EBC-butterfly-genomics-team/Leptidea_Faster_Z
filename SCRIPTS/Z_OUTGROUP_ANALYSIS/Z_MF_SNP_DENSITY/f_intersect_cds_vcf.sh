#!/bin/bash -l


module load bioinfo-tools
module load BEDTools


main_path=/home/larshook/LarsH/FastZ
snp_path=LEPTIDEA_OUTGROUPS/SNPS/MPILEUP
cds=LsinapisSweM_CDS_regions.bed
repeats=/home/larshook/LarsH/REPEATS/MASKED/P14502_103/LMC_P14502_103.FINAL-deduped-nuc.fasta.out.bed
all_windows=100kb
filter=200


for sample in P14458_103_S24 P14458_104_S25 P14458_107_S28 P14458_108_S29 P14502_103 P14502_104
do

  vcf=$main_path/$snp_path/"$sample"-filtered_"$filter".vcf.gz

  # get snps in cds, remove repeats to reduce false positives...

  zcat $vcf | grep "#" > "$vcf".temp
  
  bedtools intersect -v -a $vcf -b $repeats >> "$vcf".temp 

  bedtools intersect -c -a $main_path/$cds -b "$vcf".temp > $main_path/$snp_path/"$sample"-filtered_"$filter"_snps_per_cds.txt

  # sum in windows across genome...

  echo "chr start end snp_count" > $main_path/$snp_path/"$sample"-filtered_"$filter"_snps_per_"$all_windows"_windows.txt

  bedtools intersect -wa -wb -a $main_path/LsinapisSweM_"$all_windows".bed -b $main_path/$snp_path/"$sample"-filtered_"$filter"_snps_per_cds.txt |\
        awk '{a[$1" "$2" "$3]+=$7} END {for(i in a) print i" "a[i]}' | sort -V -k1,2 >> $main_path/$snp_path/"$sample"-filtered_"$filter"_snps_per_"$all_windows"_windows.txt

done
