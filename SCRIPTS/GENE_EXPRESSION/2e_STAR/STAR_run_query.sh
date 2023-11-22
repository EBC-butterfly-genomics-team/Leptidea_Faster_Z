#!/bin/bash -l
#SBATCH -J STAR
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 5:00:00
#SBATCH -A naiss2023-5-52


# load modules
module load bioinfo-tools 
module load star/2.7.9a
#module load cufflinks/2.2.1
module load StringTie/2.1.4
module load samtools
#module load htseq/0.6.1



READ1=${1:?msg}                         #read file names (includes path)
READ2=${2:?msg}
READ1_f=${3:?msg}						#read file names (does not include path)
READ2_f=${4:?msg}            
SRCDIR=${5:?msg}                        #read name of output folder                       


SRCDIR_INI=$(pwd)                                                                                      #remember initial path 

 
cd $SRCDIR

STAR \
      --runMode alignReads \
      --runThreadN 3 \
      --outSAMstrandField intronMotif \
      --genomeDir /home/larshook/LarsH/FastZ/RNAseq_Lsin/STAR_index/ \
      --readFilesIn $READ1 $READ2 \
      --readFilesCommand zcat \
      --outFileNamePrefix $SRCDIR/ \
      --twopassMode Basic \
      --sjdbGTFfile /proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.genes.gff \
      --sjdbGTFtagExonParentGene Parent \
      --sjdbGTFtagExonParentTranscript ID


      
# Notes (STAR, mapping step)
#--runMode alignReads				map reads
#--runThreadN N 					specifies the number of threads that will be used by the program
#--outSAMstrandField intronMotif 	adds information (to the SAM output file) required for downstream analysis with Cufflinks
#--genomeDir /path/to/index 		specifies the directory containing the pre-built genome index
#--readFilesIn /path/to/reads/sample_1.fastq /path/to/reads/sample_2.fastq 		is where you should list the FASTQ files that you wish to map
#--readFilesCommand zcat            uncompress gz files
#--outFileNamePrefix outDir 		specifies the output directory
#--twopassMode Basic				For the most sensitive novel junction discovery, run STAR in the 2-pass mode. It does not increase the number
#									of detected novel junctions, but allows to detect more splices reads mapping to novel junctions. 
#									The basic idea is to run 1st pass of STAR mapping with the usual parameters, then collect the junctions 
#									detected in the first pass, and use them as "annotated" junctions for the 2nd pass mapping.


#__________________________________________________

# Use Samtools to convert sam files into sorted bam format, and then to sort and index the bam file
cd $SRCDIR
samtools view -bS Aligned.out.sam > Aligned.out.bam
samtools sort Aligned.out.bam -o Aligned.out.sorted.bam
samtools index Aligned.out.sorted.bam


# Use Samtools to glean statistics about sorted BAM file
samtools flagstat Aligned.out.sorted.bam > STAR_flagstat.txt
samtools idxstats Aligned.out.sorted.bam > STAR_idxstats.txt


# Run StringTie

cd $SRCDIR
mkdir stringtie_2

stringtie \
      Aligned.out.sorted.bam \
      -p 3 \
      -o $SRCDIR/stringtie_2/STAR.Aligned.out.gtf \
      -A $SRCDIR/stringtie_2/STAR.Aligned.out.gene_abund.tab \
      -G /proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.genes.gff \
      -C $SRCDIR/stringtie_2/STAR.fully_covered_transcripts.gtf \
      -e 


#Notes (StringTie)
# -o 	output path/file name for the assembled transcripts GTF (default: stdout)
# -p 	number of threads (CPUs) to use (default: 1)
# -m 	minimum assembled transcript length (default: 200)
# -A <gene_abund.tab>	Gene abundances will be reported (tab delimited format) in the output file with the given name.
# -G	use reference annotation GFF to guide assembly.Output will include expressed reference transcripts.
# -C	Outputs all transcripts in the reference file that are fully covered by reads. (requires -G)
# -B	Output of Ballgown input files (*.ctab) (requires -G)
# -e	Only estimate and output assembled trancripts matching the reference transcripts given by -G (requires -G, recommended for -b/-B )

#note: runtime: 3h for largest library (3 cores); 
#      cores: only makes use of 2 or 3 cores
#      peak memory: 16GB > needs at least 3 cores


