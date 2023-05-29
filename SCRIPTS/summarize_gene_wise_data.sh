#!/bin/bash -l


# summarize all gene-wise results into one table...


module load bioinfo-tools
module load seqtk


# content...

### gene data ###
# gene_id
# gene_size
# gene position
# 4-fold GC%
# sex_bias (M, U, F)
# FPKM
# total expression level (N, L, M, H)
# female expression (N, L, M, H)
# male expression (N, L, M, H)

### chromosome data ###
# chrN: N
# Z_vs_A: Z1, Z2, Z3 or A
# anc_vs_neo: Za, Zn or A
# chr_size

### dnds data ###
# dn
# ds
# dnds
# Dn
# Ds
# Ns
# Ss

### polymorphism data ###
# Pn
# Ps

### pi 0/4-fold ###
# pi0
# pi4

## Tajima's D
# TajD



######

main_path=/home/larshook/LarsH/FastZ
msa_path=DN_DS/FILTERED_ALIGNMENTS

lep_gff=/proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.genes.gff
assembly=/proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/LsinapisSweM_chr.fasta
gc_perc=/home/larshook/LarsH/FastZ/PI/LsinapisSweM_4_fold_gc.txt
sbg=RNAseq_Lsin/DESEQ2/LFC_Adult.txt

exp_class=/home/larshook/SCRIPTS/FastZ/classified_expressed_genes.txt
exp_levels=/home/larshook/SCRIPTS/FastZ/mean_FPKM_expressed_genes.txt
circos_path=/home/larshook/SCRIPTS/CIRCOS/bm_ls



######

# extact some data from result files...

# make a list of all gene numbers in the analysis...
sed 's/gene_//;s/\.best.filtered.phy//' "$main_path"/"$msa_path"/gene_list.txt > "$main_path"/analyzed_genes_list.txt

# make file with chr size info...
seqtk comp "$assembly" | awk '{print $1, $2}' > lep_sin_chr_sizes.txt

# make file with dnds, dn, ds, Dn and Ds from paml output...
printf "gene_id\tt\tN\tS\tdN/dS\tdN\tdS\tDn\tDs\n" > "$main_path"/all_paml_output.txt
cat "$main_path"/DN_DS/genes_*_block.out | grep -A 3 "#1: LSSWET\|branch" |\
	awk '{if ($1=="#1:" || $1=="7..1") print $0}' | paste - - | awk -v OFS="\t" '{print $2, $4, $5, $6, $7, $8, $9, $10, $11}' >> "$main_path"/all_paml_output.txt

# get coordinates for ancestral and neo Z1 from collinearity...
neo_end=$(grep -w "ls1" "$circos_path"/bm_ls.collinearity.links | grep -w "bm17" | awk '{print ($5>$6?$5:$6)}' | sort -Vr | awk 'NR==1 {print $0}')
anc_start=$(grep -w "ls1" "$circos_path"/bm_ls.collinearity.links | grep -w "bm1" | awk '{print ($5<$6?$5:$6)}' | sort -V | awk 'NR==1 {print $0}')

neo_anc_break=$(awk -v neo_end="$neo_end" -v anc_start="$anc_start" 'BEGIN{OFMT="%.0f"; print (neo_end+anc_start)/2}')



######


# make summary table file...

printf "gene_id\tgene_size\tgene_position\tgc4fold\tsex_bias\tFPKM\ttotal_exp_level\tfemale_exp_level\tmale_exp_level\tchrN\tZ_vs_A\tanc_vs_neo\tchr_size\tdn\tds\tdnds\tDn\tDs\tNs\tSs\tPn\tPs\tpi0\tpi4\tTajD\n" > "$main_path"/lep_sin_dnds_data.txt

# loop through analyzed_genes_list.txt and create variables for each column for each gene and print to lep_sin_dnds_data.txt...

while read line
do

  gene_id=$(printf "LSSWET$line")
  gene_size=$(grep -w "LSSWEG$line" $lep_gff | awk '{if ($3=="gene") print $5-$4}')
  gene_position=$(grep -w "LSSWEG$line" $lep_gff | awk '{if ($3=="gene") printf("%.0f\n", ($4+$5)/2)}')
  gc4fold=$(grep -w "LSSWET$line" $gc_perc | awk '{print $2}') 
  sex_bias=$(grep -w "LSSWEG$line" "$main_path"/"$sbg" | awk '{if ($3>1 && $7<0.05) print "MBG"; else if ($3<-1 && $7<0.05) print "FBG"; else print "UBG"}')
  FPKM=$(grep -w "LSSWEG$line" $exp_levels | awk '{print $2}')
  total_exp_level=$(grep -w "LSSWEG$line" $exp_class | grep -w "total" | awk '{print $2}')
  female_exp_level=$(grep -w "LSSWEG$line" $exp_class | grep -w "female" | awk '{print $2}')
  male_exp_level=$(grep -w "LSSWEG$line" $exp_class | grep -w "male" | awk '{print $2}')

  chrN=$(grep -w "LSSWEG$line" $lep_gff | awk '{if ($3=="gene") print $1}')
  Z_vs_A=$(grep -w "LSSWEG$line" $lep_gff | awk '{if ($1=="HiC_scaffold_1" || $1=="HiC_scaffold_6" || $1=="HiC_scaffold_18") print "Z"; else print "A"}' | uniq)
  anc_vs_neo=$(grep -w "LSSWEG$line" $lep_gff |\
	awk -v neo_anc_break=$neo_anc_break '{if (($1=="HiC_scaffold_1" && $5<neo_anc_break) || $1=="HiC_scaffold_6" || $1=="HiC_scaffold_18") print "neo Z";
		else if ($1=="HiC_scaffold_1" && $4>neo_anc_break) print "anc Z"; else print "A"}' | uniq)
  chr_size=$(grep -w "$chrN" lep_sin_chr_sizes.txt | awk '{print $2}')

  dn=$(grep -w "LSSWET$line" "$main_path"/all_paml_output.txt | awk '{print $6}')
  ds=$(grep -w "LSSWET$line" "$main_path"/all_paml_output.txt | awk '{print $7}')
  dnds=$(grep -w "LSSWET$line" "$main_path"/all_paml_output.txt | awk '{print $5}')
  Dn=$(grep -w "LSSWET$line" "$main_path"/all_paml_output.txt | awk '{print $8}')
  Ds=$(grep -w "LSSWET$line" "$main_path"/all_paml_output.txt | awk '{print $9}')
  Ns=$(grep -w "LSSWET$line" "$main_path"/all_paml_output.txt | awk '{print $3}')
  Ss=$(grep -w "LSSWET$line" "$main_path"/all_paml_output.txt | awk '{print $4}')

  Pn=$(grep -w "LSSWET$line" "$main_path"/PN_PS/mk.tsv | awk '{print $2}')
  Ps=$(grep -w "LSSWET$line" "$main_path"/PN_PS/mk.tsv | awk '{print $3}')

  pi0=$(grep -w "LSSWEG$line" "$main_path"/PI/gene_wise_0-fold_pi.txt | awk '{print $2}')
  pi4=$(grep -w "LSSWEG$line" "$main_path"/PI/gene_wise_4-fold_pi.txt | awk '{print $2}')

  TajD=$(grep -w "LSSWEG$line" "$main_path"/TAJIMAS_D/gene_wise_Tajimas_D.txt | awk '{print $2}')

  printf "$gene_id\t$gene_size\t$gene_position\t$gc4fold\t$sex_bias\t$FPKM\t$total_exp_level\t$female_exp_level\t$male_exp_level\t$chrN\t$Z_vs_A\t$anc_vs_neo\t$chr_size\t$dn\t$ds\t$dnds\t$Dn\t$Ds\t$Ns\t$Ss\t$Pn\t$Ps\t$pi0\t$pi4\t$TajD\n" >> /home/larshook/LarsH/FastZ/lep_sin_dnds_data.txt

done < "$main_path"/analyzed_genes_list.txt


# rename scaffolds to chr and Z for plotting...

sed -i 's/HiC_scaffold_/Chr /' /home/larshook/LarsH/FastZ/lep_sin_dnds_data.txt
sed -i 's/Chr 1\t/Z1\t/;s/Chr 6/Z2/;s/Chr 18/Z3/' /home/larshook/LarsH/FastZ/lep_sin_dnds_data.txt
