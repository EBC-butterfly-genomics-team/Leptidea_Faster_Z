#!/bin/bash -l
#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH -J make_variant_fasta

# snp mask and make alternative snp assembly fasta based on vcf

module load bioinfo-tools
module load GATK/4.1.1.0
module load java/sun_jdk1.8.0_151
module load samtools

vcf=LsinapisSweM-Swe-sin_final_all-filtered_q30_snps.vcf
reference=LsinapisSweM
main_path=/home/larshook/LarsH/FastZ
vcf_path=Snps/HaplotypeCaller/round2


cp $main_path/$vcf_path/$vcf $SNIC_TMP
cp /proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/$reference.fasta $SNIC_TMP

cd $SNIC_TMP

samtools faidx $reference.fasta
gatk --java-options "-Xmx6g" CreateSequenceDictionary -R $reference.fasta
gatk --java-options "-Xmx6g" IndexFeatureFile -F "$vcf"

# make snp masked fasta...

gatk --java-options "-Xmx6g" FastaAlternateReferenceMaker \
	-R $reference.fasta \
	-O "$reference"_snp_masked.fasta \
	-V $vcf \
	--snp-mask $vcf \
	--snp-mask-priority

# reformat header to original format...

sed 's/:/ /' "$reference"_snp_masked.fasta | awk '{if (NF>1) print ">"$2; else print}' > "$reference"_snp_masked.temp
mv "$reference"_snp_masked.temp "$reference"_snp_masked.fasta

# make variant fasta...

gatk --java-options "-Xmx6g" FastaAlternateReferenceMaker \
        -R $reference.fasta \
        -O "$reference"_variant.fasta \
        -V $vcf

# reformat header to original format...

sed 's/:/ /' "$reference"_variant.fasta | awk '{if (NF>1) print ">"$2; else print}' > "$reference"_variant.temp
mv "$reference"_variant.temp "$reference"_variant.fasta



cp "$reference"_snp_masked.fasta "$reference"_variant.fasta $main_path
