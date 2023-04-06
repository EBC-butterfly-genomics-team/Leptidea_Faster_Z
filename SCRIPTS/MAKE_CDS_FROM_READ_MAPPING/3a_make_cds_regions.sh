#!/bin/bash


# create region file of CDS from gff for use with samtools consensus

out_dir=/home/larshook/LarsH/FastZ



# Swe sinapis male...

awk '{if ($3=="CDS") print $1":"$4"-"$5}' /proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.genes.gff |\
	sort -V > "$out_dir"/LsinapisSweM_CDS_regions.txt

sed 's/ID=//;s/-RA/ /' /proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.genes.gff |\
	awk '{if ($3=="CDS") print $1":"$4"-"$5, $7, $9}' |\
		sort -V > "$out_dir"/LsinapisSweM_CDS_regions_and_gene_info.txt



# Swe sinapis female...

awk '{if ($3=="CDS") print $1":"$4"-"$5}' /home/larshook/LarsH/MAKER/SIN_Swe_female/round4/LMC_P14502_104.FINAL-deduped-nuc.all.maker.genes.domain.putative_function.gff |\
        sort -V > "$out_dir"/LsinapisSweF_CDS_regions.txt

sed 's/ID=//;s/-RA/ /' /home/larshook/LarsH/MAKER/SIN_Swe_female/round4/LMC_P14502_104.FINAL-deduped-nuc.all.maker.genes.domain.putative_function.gff |\
        awk '{if ($3=="CDS") print $1":"$4"-"$5, $7, $9}' |\
                sort -V > "$out_dir"/LsinapisSweF_CDS_regions_and_gene_info.txt

