#!/bin/bash


# create region file of CDS from gff for use with samtools consensus

cd /home/larshook/LarsH/FastZ/DN_DS


awk '{if ($3=="CDS") print $1":"$4"-"$5}' /proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.genes.gff |\
	sort -V > LsinapisSweM_CDS_regions.txt

sed 's/ID=//;s/-RA/ /' /proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.genes.gff |\
	awk '{if ($3=="CDS") print $1":"$4"-"$5, $7, $9}' |\
		sort -V > LsinapisSweM_CDS_regions_and_gene_info.txt
