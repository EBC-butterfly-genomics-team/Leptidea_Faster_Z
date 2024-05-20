#!/bin/bash -l


# summarize all gene-wise results into one table...


module load bioinfo-tools
module load seqtk


# content...

### gene data ###
# gene_id
# gene_size
# gene position

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

######

main_path=/home/larshook/LarsH/FastZ
msa_path=DN_DS/FILTERED_ALIGNMENTS

lep_gff=/proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.genes.gff
assembly=/proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/LsinapisSweM_chr.fasta

circos_path=/home/larshook/SCRIPTS/CIRCOS/bm_ls



######

# extact some data from result files...

# make a list of all gene numbers in the analysis...
sed 's/gene_//;s/\.best.filtered.phy//' "$main_path"/"$msa_path"/gene_list.txt > "$main_path"/analyzed_genes_list.txt

# make file with chr size info...
seqtk comp "$assembly" | awk '{print $1, $2}' > lep_sin_chr_sizes.txt



# make file with dnds, dn and ds from paml output...
# 7..1 = Lsin
# 9..3 = Lmor
# 6..5 = Ldup

printf "gene_id\tLsin_dndn\tLsin_dn\tLsin_ds\tLmor_dndn\tLmor_dn\tLmor_ds\tLdup_dndn\tLdup_dn\tLdup_ds\n" > "$main_path"/mb_paml_output.txt
cat "$main_path"/DN_DS/MULTI_BRANCH/genes_*_block.out | grep -A 9 "#1: LSSWET\|branch" |\
	awk '{if ($1=="#1:" || $1=="7..1" || $1=="9..3" || $1=="6..5") print $0}' | paste - - - - |\
		awk -v OFS="\t" '{print $2, $7, $8, $9, $16, $17, $18, $25, $26, $27}' >> "$main_path"/mb_paml_output.txt



# get coordinates for ancestral and neo Z1 from collinearity...
neo_end=$(grep -w "ls1" "$circos_path"/bm_ls.collinearity.links | grep -w "bm17" | awk '{print ($5>$6?$5:$6)}' | sort -Vr | awk 'NR==1 {print $0}')
anc_start=$(grep -w "ls1" "$circos_path"/bm_ls.collinearity.links | grep -w "bm1" | awk '{print ($5<$6?$5:$6)}' | sort -V | awk 'NR==1 {print $0}')

neo_anc_break=$(awk -v neo_end="$neo_end" -v anc_start="$anc_start" 'BEGIN{OFMT="%.0f"; print (neo_end+anc_start)/2}')



######


# make summary table file...

printf "gene_id\tgene_size\tgene_position\tchrN\tZ_vs_A\tanc_vs_neo\tchr_size\tLsin_dnds\tLsin_dn\tLsin_ds\tLmor_dnds\tLmor_dn\tLmor_ds\tLdup_dnds\tLdup_dn\tLdup_ds\n" > "$main_path"/lep_sin_mb_dnds_data.txt

# loop through analyzed_genes_list.txt and create variables for each column for each gene and print to lep_sin_dnds_data.txt...

while read line
do

  gene_id=$(printf "LSSWET$line")
  gene_size=$(grep -w "LSSWEG$line" $lep_gff | awk '{if ($3=="gene") print $5-$4}')
  gene_position=$(grep -w "LSSWEG$line" $lep_gff | awk '{if ($3=="gene") printf("%.0f\n", ($4+$5)/2)}')

  chrN=$(grep -w "LSSWEG$line" $lep_gff | awk '{if ($3=="gene") print $1}')
  Z_vs_A=$(grep -w "LSSWEG$line" $lep_gff | awk '{if ($1=="HiC_scaffold_1" || $1=="HiC_scaffold_6" || $1=="HiC_scaffold_18") print "Z"; else print "A"}' | uniq)
  anc_vs_neo=$(grep -w "LSSWEG$line" $lep_gff |\
	awk -v neo_anc_break=$neo_anc_break '{if (($1=="HiC_scaffold_1" && $5<neo_anc_break) || $1=="HiC_scaffold_6" || $1=="HiC_scaffold_18") print "neo Z";
		else if ($1=="HiC_scaffold_1" && $4>neo_anc_break) print "anc Z"; else print "A"}' | uniq)
  chr_size=$(grep -w "$chrN" lep_sin_chr_sizes.txt | awk '{print $2}')
  Lsin_dnds=$(grep -w "LSSWET$line" "$main_path"/mb_paml_output.txt | awk '{print $2}')
  Lsin_dn=$(grep -w "LSSWET$line" "$main_path"/mb_paml_output.txt | awk '{print $3}')
  Lsin_ds=$(grep -w "LSSWET$line" "$main_path"/mb_paml_output.txt | awk '{print $4}')
  Lmor_dnds=$(grep -w "LSSWET$line" "$main_path"/mb_paml_output.txt | awk '{print $5}')
  Lmor_dn=$(grep -w "LSSWET$line" "$main_path"/mb_paml_output.txt | awk '{print $6}')
  Lmor_ds=$(grep -w "LSSWET$line" "$main_path"/mb_paml_output.txt | awk '{print $7}')
  Ldup_dnds=$(grep -w "LSSWET$line" "$main_path"/mb_paml_output.txt | awk '{print $8}')
  Ldup_dn=$(grep -w "LSSWET$line" "$main_path"/mb_paml_output.txt | awk '{print $9}')
  Ldup_ds=$(grep -w "LSSWET$line" "$main_path"/mb_paml_output.txt | awk '{print $10}')

  printf "$gene_id\t$gene_size\t$gene_position\t$chrN\t$Z_vs_A\t$anc_vs_neo\t$chr_size\t$Lsin_dnds\t$Lsin_dn\t$Lsin_ds\t$Lmor_dnds\t$Lmor_dn\t$Lmor_ds\t$Ldup_dnds\t$Ldup_dn\t$Ldup_ds\n" >> /home/larshook/LarsH/FastZ/lep_sin_mb_dnds_data.txt

done < "$main_path"/analyzed_genes_list.txt


# rename scaffolds to chr and Z for plotting...

sed -i 's/HiC_scaffold_/Chr /' /home/larshook/LarsH/FastZ/lep_sin_mb_dnds_data.txt
sed -i 's/Chr 1\t/Z1\t/;s/Chr 6/Z2/;s/Chr 18/Z3/' /home/larshook/LarsH/FastZ/lep_sin_mb_dnds_data.txt
