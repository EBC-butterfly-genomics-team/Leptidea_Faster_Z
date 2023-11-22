#!/bin/bash -l
#SBATCH -J condetri
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -A snic2022-5-34

## --mail-user lars.hook@ebc.uu.se
## --mail-type=ALL

# load modules
module load bioinfo-tools
module load FastQC/0.11.8


READ1=${1:?msg}                      #read name of files containing paired reads (with path)
READ2=${2:?msg}
REFNAME1=${3:?msg}                   #read name of files containing paired reads (no path)
REFNAME2=${4:?msg}
SRCDIR=${5:?msg}                     #output folder location

cd $SRCDIR

#unzip files (or just copy to scratch if already unzipped)

#READ1_fq=${READ1#/proj/b2014034/nobackup/Lars/RNAtesting/Drosophila/1d_QC_filterPolyA_prinseq/polyA_OUT-masked-}
#READ2_fq=${READ2#/proj/b2014034/nobackup/Lars/RNAtesting/Drosophila/1d_QC_filterPolyA_prinseq/polyA_OUT-masked-}

REFNAME1_fq=${REFNAME1#polyA_OUT-masked-}
REFNAME2_fq=${REFNAME2#polyA_OUT-masked-}

cp $READ1 $SNIC_TMP/${REFNAME1_fq}
cp $READ2 $SNIC_TMP/${REFNAME2_fq}



cd $SNIC_TMP                                           # use local scratch disk


# Usage:   perl condetri.pl -fastq1=file1 [-fastq2=file2 -prefix=s -cutfirst=i -cutlast=i -rmN -notrim -hq=i -lq=i -frac=[0,1] -lfrac=[0,1] -minlen=i -mh=i -ml=i -sc=i -pb=s]
PRXOUT=condetri-${REFNAME1_fq%_1.fastq}              #output file(s) prefix

perl /home/larshook/Lars/RNAseq_Lsin_2020/SCRIPTS/1e_QC_condetri/condetri.pl \
            -fastq1=$SNIC_TMP/${REFNAME1_fq} \
            -fastq2=$SNIC_TMP/${REFNAME2_fq} \
            -prefix=$PRXOUT \
            -cutfirst=0 \
            -cutlast=0 \
            -rmN \
            -hq=30 \
            -lq=0 \
            -frac=0.8 \
            -minlen=30 \
            -mh=1 \
            -ml=1 \
            -sc=33
                        
fastqc -o $SRCDIR/ ${PRXOUT}_trim1.fastq        #run FastQC on filtered fastq files
fastqc -o $SRCDIR/ ${PRXOUT}_trim2.fastq

              

# notes:  
# -fastq1		first file of paired-end reads 
# -fastq2		second file of paired-end reads 
# -prefix		string Prefix for the output file(s) 
# -cutfirst=i 	Remove i bases from 5'-end before any trimming 
# -cutlast=i 	Remove i bases from the 3'-end before any trimming  
# -rmN 			Remove Ns from 5' end before any trimming 
# -hq:			Bases are removed from the 3'-end if the quality score is lower  
#				than some hq quality threshold;	 
# -ml			Max number of lq bases allowed within a stretch of hq bases 
# -lq	 		Low quality threshold (read filtered out if there is a base in the
#				read with quality<lq)
# -mh 			When i consecutive hq bases is reached, the trimming stops 
# -frac=[0,1]   Fraction of read that must exceed hq after quality trimming 
# -lfrac=[0,1]  Maximum fraction of bases with qual<lq after quality trimming[0]
# -minlen       Min allowed read length
# -sc=i			Illumina scoring table, Score=ASCII-sc. Used to be 64 for Illumina/Solexa. 
# 				33 is Sanger standard, and also used for newer Illumina data (1.8+)

multiqc $SRCDIR/                  # creates single report
multiqc $SNIC_TMP                  # creates single report



cp -r multiqc_data $SRCDIR/       #copy results from local disk back to home directory
cp * $SRCDIR/
              



#note: runtime when using project folder (for one 3G fastq file): 1h
#                                                    largest file: 8h


