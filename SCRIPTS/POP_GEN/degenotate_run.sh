#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J degenotate

module load bioinfo-tools
module load biopython
module load BEDTools


main_path=/home/larshook/LarsH/FastZ

# get 4- and 0-fold degenerate sites...

python /home/larshook/LarsH/SOFTWARE/Degenotate/degenotate.py \
	-a "$main_path"/DN_DS/LsinapisSweM.all.maker.genes_rebuilt.gff \
	-g /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/LsinapisSweM.fasta \
	-o "$main_path"/PI \
	--overwrite


# get chromosome sites only...

for i in {1..29}
do 
  grep -w HiC_scaffold_"$i" "$main_path"/PI/degeneracy-all-sites.bed
done > "$main_path"/PI/degeneracy-all-sites_chr.bed


# restrict degenotate output to updated CDS coords...

awk -v OFS="\t" '{if ($3=="CDS") print $1, $4-1, $5}' "$main_path"/DN_DS/recalculated_LsinapisSweM.all.maker.genes.gff |\
	bedtools intersect -a "$main_path"/PI/degeneracy-all-sites_chr.bed -b - > "$main_path"/PI/degeneracy-all-sites_CDS.bed


# make 4- and 0-fold sites files for pixy...

awk -v OFS="\t" '{if ($5=="0") print $1, $3}' "$main_path"/PI/degeneracy-all-sites_CDS.bed > "$main_path"/PI/LsinapisSweM_0-fold_sites.txt
awk -v OFS="\t" '{if ($5=="4") print $1, $3}' "$main_path"/PI/degeneracy-all-sites_CDS.bed > "$main_path"/PI/LsinapisSweM_4-fold_sites.txt

# make all CDS sites file for pixy...

awk -v OFS="\t" '{print $1, $3}' "$main_path"/PI/degeneracy-all-sites_CDS.bed > "$main_path"/PI/LsinapisSweM_N-fold_sites.txt
