#!/bin/bash -l
#SBATCH -J fq_masker
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 05:00:00
#SBATCH -A snic2020-5-20
#SBATCH --mail-user lars.hook@ebc.uu.se
#SBATCH --mail-type=ALL

# load modules
module load bioinfo-tools
module load Fastx/0.0.14
module load FastQC/0.11.8



READ1=${1:?msg}                      #read name of files containing paired reads (with path)
READ2=${2:?msg}
REFNAME1=${3:?msg}                   #read name of files containing paired reads (no path)
REFNAME2=${4:?msg}
SRCDIR=${5:?msg}                             #location output folder


#unzip fq.gz files
gunzip -c $READ1>$SNIC_TMP/READ1.fastq
gunzip -c $READ2>$SNIC_TMP/READ2.fastq

#cd $SNIC_TMP                                # use local scratch disk



# fastq_masker is being used to mask low quality nucleotides

REF_PREF=masked

 
# mask nucleotides with quality below 20
# each file (paired end) masked independently

fastq_masker \
        -v \
        -q 20 \
        -r N \
        -z \
        -i $SNIC_TMP/READ1.fastq \
        -o $SRCDIR/${REF_PREF}-${REFNAME1}


fastq_masker \
        -v \
        -q 20 \
        -r N \
        -z \
        -i $SNIC_TMP/READ2.fastq \
        -o $SRCDIR/${REF_PREF}-${REFNAME2}


# notes:
#   [-q N]       = Quality threshold - nucleotides with lower quality will be masked
#                  Default is 10.
#   [-r C]       = Replace low-quality nucleotides with character C. Default is 'N'
#   [-z]         = Compress output with GZIP.
#   [-i INFILE]  = FASTQ input file. default is STDIN.
#   [-o OUTFILE] = FASTQ output file. default is STDOUT.
#   [-v]         = Verbose - report number of sequences.
#                  If [-o] is specified,  report will be printed to STDOUT.
#                  If [-o] is not specified (and output goes to STDOUT),
#                  report will be printed to STDERR.




#perform FastQC analysis
fastqc -o $SRCDIR/ $SRCDIR/${REF_PREF}-${REFNAME1}
fastqc -o $SRCDIR/ $SRCDIR/${REF_PREF}-${REFNAME2}



#multiqc $SRCDIR/                  # creates single report


#note: runtime (for 2x fq.gz files; 1.9-3.1 Gb each): up to 1h40m

