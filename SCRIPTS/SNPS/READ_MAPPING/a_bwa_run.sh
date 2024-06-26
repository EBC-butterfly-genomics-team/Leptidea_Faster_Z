#!/bin/bash

for sample in \
	Swe-sin-1C_GCCAAT \
	Swe-sin-2C_CAGATC \
	Swe-sin-31C_ACTTGA \
	Swe-sin-32C_GATCAG \
	Swe-sin-61C_TAGCTT \
	Swe-sin-101C_AGTTCC \
	Swe-sin-102C_ATGTCA \
	Swe-sin-62C_GGCTAC \
	Swe-sin-91C_CTTGTA \
	Swe-sin-92C_AGTCAA \
	Ire-juv-1C_CCGTCC \
	Ire-juv-21C_GTGAAA \
	Ire-juv-22C_GTGGCC \
	Ire-juv-2C_GTCCGC \
	Ire-juv-41C_GTTTCG \
	Ire-juv-42C_CGTACG \
	Ire-juv-61C_GAGTGG \
	Ire-juv-62C_ACTGAT \
	Ire-juv-81C_ATTCCT \
	Ire-juv-82C_ATCACG
do
  sbatch bwa_query.sh $sample
  sleep 2
done

cd /home/larshook/LarsH/FastZ/

mkdir -p Snps
cd Snps
mkdir -p Read_mapping
