#!/bin/bash


# update gff with recalculated coordinates...


main_path=/home/larshook/LarsH/FastZ/DN_DS



# genes that for some reason don't have complete orfs, from recalculate_cds_coordinates.sh error log...

error_log=...



# make list with all genes to extract from original gff...

grep Error "$error_log" |\
	awk '{print $6}' |\
		sed "s/'ORFS\/LSREFT//;s/.fasta'//" |\
			awk '{print "LSSWET"$1}' |\
				uniq > "$main_path"/coords_to_update.txt

awk 'gene=NF {print $gene}' "$main_path"/recalculated_gene_coordinates.txt |\
	sed 's/ID=//' |\
		uniq >> "$main_path"/coords_to_update.txt



# use list to make updated gff without recalculated coordinates...

grep -vwf "$main_path"/coords_to_update.txt /proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.genes.gff > "$main_path"/recalculated_LsinapisSweM.all.maker.genes.gff



# add recalculated coords to new gff, skip non-orf...

cat "$main_path"/recalculated_gene_coordinates.txt |\
	awk '{if ($4!=".") print}' >> "$main_path"/recalculated_LsinapisSweM.all.maker.genes.gff



# restructure gff so it works with Degenotate... (gene structure is not "correct" in Maker one-transcript gff)
# remove transcript name from CDS row, this should only be on mRNA rows not on both, structure should be: CDS -> parent -> mRNA -> parent -> gene
# also remove ";" linebreak

sed 's/;/ /g' "$main_path"/recalculated_LsinapisSweM.all.maker.genes.gff |\
	awk -v OFS="\t" '{if ($3=="CDS") print $1, $2, $3, $4, $5, $6, $7, $8, $10; else print}' |\
		sed 's/ /;/' |\
			awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' > "$main_path"/LsinapisSweM.all.maker.genes_rebuilt.gff


# make 1-based gene bed for pixy, gff is 1-based...
# "Note that pixy expects intervals to be one indexed, not zero indexed as many typical BED files are.
# If your BED file is zero indexed, you’ll need to convert the intervals (e.g. by adding 1 to the start and end)."


for i in {1..29}; do grep -w HiC_scaffold_"$i" LsinapisSweM.all.maker.genes_rebuilt.gff; done |\
	sed 's/ID=//;s/;/ /' |\
		awk -v OFS="\t" '{if ($3=="gene") print $1, $4, $5, $9}' > /home/larshook/LarsH/FastZ/PI/LsinapisSweM.all.maker.genes.bed

