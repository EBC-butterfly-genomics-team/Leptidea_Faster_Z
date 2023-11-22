#!/bin/bash -l
#SBATCH -J PolyA_prinseq
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 20:00:00
#SBATCH -A snic2020-5-20
#SBATCH --mail-user lars.hook@ebc.uu.se
#SBATCH --mail-type=ALL

# load modules
module load bioinfo-tools
module load prinseq/0.20.4
module load FastQC/0.11.8
#module load MultiQC/0.9



READ1=${1:?msg}                      #read name of files containing paired reads (with path)
READ2=${2:?msg}
REFNAME1=${3:?msg}                   #read name of files containing paired reads (no path)
REFNAME2=${4:?msg}
SRCDIR=${5:?msg}                     #location output folder


#unzip fq.gz files
gunzip -c $READ1>$SNIC_TMP/READ1.fastq
gunzip -c $READ2>$SNIC_TMP/READ2.fastq

#cd $SNIC_TMP                                # use local scratch disk



# prinseq is being used only to trim poly(A) tail, Poly(N) and/or low complexity segments

REF_PREF=polyA_OUT
REF_PREF2=prinseq_REJECTED


# Trim Poly(A), Poly(N), low complexity segments, and high N content reads
#       - prinseq trimming is done in the following order:
#       :trim ends with qual<30
#		:trim poly(A) on both ends
#       :trim poly(N) on both ends
#       :filter low complexity segments using DUST method
#       	note: '-lc_method dust' with '-lc_threshold 7' option filters T{19} and longer
#   		(The DUST approach is adapted from the algorithm used to mask low-complexity regions 
#   		during BLAST search preprocessing)
#       :minimum length: 30
#       :Filter reads with high N content(>10%)


prinseq \
           -verbose \
           -fastq $SNIC_TMP/READ1.fastq \
           -fastq2 $SNIC_TMP/READ2.fastq \
           -out_good $SRCDIR/${REF_PREF}-${REFNAME1%_R1_001_val_1.fq.gz} \
           -out_bad $SRCDIR/${REF_PREF2}-${REFNAME1%_R1_001_val_1.fq.gz} \
           -log $SRCDIR/${REFNAME1%_R1_001_val_1.fq.gz}-prinseq.log \
           -graph_data $SRCDIR/${REFNAME1%_R1_001_val_1.fq.gz}-prinseq.gd \
           -graph_stats ld,gc,qd, ns, pt, ts, aq, de, sc, dn \
           -qual_noscale \
           -exact_only \
           -min_len 30 \
           -ns_max_p 10 \
           -trim_tail_left 2 \
           -trim_tail_right 2 \
           -lc_method dust \
           -lc_threshold 7 \
           -trim_ns_left 1 \
           -trim_ns_right 1 \
           -trim_qual_left 30 \
           -trim_qual_right 30

         

#perform analysis after trimming Poly(A) etc

#prinseq \
#           -verbose \
#           -fastq $SRCDIR/${REF_PREF}-${REFNAME1%_1_val_1.fq.gz}_1.fastq \
#           -fastq2 $SRCDIR/${REF_PREF}-${REFNAME2%_2_val_2.fq.gz}_2.fastq \
#           -out_good null \
#           -out_bad null \
#           -log $SRCDIR/${REFNAME1%_1_val_1.fq.gz}-prinseq-FINAL.log \
#           -graph_data $SRCDIR/${REFNAME1%_1_val_1.fq.gz}-prinseq-FINAL.gd \
#           -graph_stats ld,gc,qd, ns, pt, ts, aq, de, sc, dn \
#           -qual_noscale \
#           -exact_only
       


#perform FastQC analysis
fastqc -o $SRCDIR/ $SRCDIR/${REF_PREF}-${REFNAME1%_R1_001_val_1.fq.gz}_1.fastq
fastqc -o $SRCDIR/ $SRCDIR/${REF_PREF}-${REFNAME2%_R2_001_val_2.fq.gz}_2.fastq

#multiqc $SRCDIR/                  # creates single report




           
# notes (prinseq):
#-verbose
#            Prints status and info messages during processing.
#-fastq <file>
#            Input file in FASTQ format that contains the sequence and
#            quality data. Use stdin instead of a file name to read from
#            STDIN (-fasta stdin). This can be useful to process compressed
#            files using Unix pipes.
#-fastq2 <file>
#            For paired-end data only. Input file in FASTQ format that
#            contains the sequence and quality data. The sequence identifiers
#            for two matching paired-end sequences in separate files can be
#            marked by /1 and /2, or _L and _R, or _left and _right, or must
#            have the exact same identifier in both input files. The input
#            sequences must be sorted by their sequence identifiers.
#            Singletons are allowed in the input files.
#-out_good <string>
#            By default, the output files are created in the same directory
#            as the input file containing the sequence data with an
#            additional "_prinseq_good_XXXX" in their name (where XXXX is
#            replaced by random characters to prevent overwriting previous
#            files). To change the output filename and location, specify the
#            filename using this option. The file extension will be added
#            automatically (either .fasta, .qual, or .fastq). For paired-end
#            data, filenames contain additionally "_1", "_1_singletons",
#            "_2", and "_2_singletons" before the file extension. Use
#            "-out_good null" to prevent the program from generating the
#            output file(s) for data passing all filters. Use "-out_good
#            stdout" to write data passing all filters to STDOUT (only for
#            FASTA or FASTQ output files).
#            Example: use "file_passed" to generate the output file
#            file_passed.fasta in the current directory
#-out_bad <string>
#            By default, the output files are created in the same directory
#            as the input file containing the sequence data with an
#            additional "_prinseq_bad_XXXX" in their name (where XXXX is
#            replaced by random characters to prevent overwriting previous
#            files). To change the output filename and location, specify the
#            filename using this option. The file extension will be added
#            automatically (either .fasta, .qual, or .fastq). For paired-end
#            data, filenames contain additionally "_1" and "_2" before the
#            file extension. Use "-out_bad null" to prevent the program from
#            generating the output file(s) for data not passing any filter.
#            Use "-out_bad stdout" to write data not passing any filter to
#            STDOUT (only for FASTA or FASTQ output files).
#            Example: use "file_filtered" to generate the output file
#            file_filtered.fasta in the current directory
#            Example: "-out_good stdout -out_bad null" will write data
#            passing filters to STDOUT and data not passing any filter will
#            be ignored.
#-log <file>
#            Log file to keep track of parameters, errors, etc. The log file
#            name is optional. If no file name is given, the log file name
#            will be "inputname.log". If the log file already exists, new
#            content will be added to the file.
#-graph_data <file>
#            File that contains the necessary information to generate the
#            graphs similar to the ones in the web version. The file name is
#            optional. If no file name is given, the file name will be
#            "inputname.gd". If the file already exists, new content will
#            overwrite the file. Use "-out_good null -out_bad null" to
#            prevent generating any additional outputs. (See below for more
#            options related to the graph data.)
#            The graph data can be used as input for the prinseq-graphs.pl
#            file to generate the PNG graph files or an HTML report file. If
#            you have trouble installing the required prinseq-graphs.pl
#            modules or want to see an output example report, upload the
#            graph data file at: http://edwards.sdsu.edu/prinseq/ -> Choose
#            "Get Report"
#-graph_stats <string>
#            Use this option to select what statistics should be calculated
#            and included in the graph_data file. This is useful if you e.g.
#            do not need sequence complexity information, which requires a
#            lot of computation. Requires to have graph_data specified.
#            Default is all selected.
#
#            Allowed option are (separate multiple by comma with no spaces):
#            ld (Length distribution), gc (GC content distribution), qd (Base
#            quality distribution), ns (Occurence of N), pt (Poly-A/T tails),
#            ts (Tag sequence check), aq (Assembly quality measure), de
#            (Sequence duplication - exact only), da (Sequence duplication -
#            exact + 5'/3'), sc (Sequence complexity), dn (Dinucleotide odds
#            ratios, includes the PCA plots)
#            Example use: -graph_stats ld,gc,qd,de
#-qual_noscale
#            Use this option if all your sequences are shorter than 100bp as
#            they do not require to scale quality data to 100 data points in
#            the graph. By default, quality scores of sequences shorter than
#            100bp or longer than 100bp are fit to 100 data points. (To
#            retrieve this information and calculate the graph data would
#            otherwise require to parse the data two times or store all the
#            quality data in memory.)
#-exact_only
#            Use this option to check for exact (forward and reverse)
#            duplicates only when generating the graph data. This allows to
#            keep the memory requirements low for large input files and is
#            faster. This option will automatically be applied when using
#            -derep options 1 and/or 4 only. Specify option -derep 1 or
#            -derep 4 if you do not want to apply both at the same time.
#-min_len <integer>
#            Filter sequence shorter than min_len.
#-ns_max_p <integer>
#            Filter sequence with more than ns_max_p percentage of Ns.
#-lc_method <string>
#            Method to filter low complexity sequences. The current options
#            are "dust" and "entropy". Use "-lc_method dust" to calculate the
#            complexity using the dust method.
#            note: '-lc_method dust' with '-lc_threshold 7' option filters T{19} and longer
#-lc_threshold <integer>
#            The threshold value (between 0 and 100) used to filter sequences
#            by sequence complexity. The dust method uses this as maximum
#            allowed score and the entropy method as minimum allowed value.
#-trim_tail_right <integer>
#            Trim poly-A/T tail with a minimum length of trim_tail_right at
#            the 3'-end.
#-trim_ns_left <integer>
#            Trim poly-N tail with a minimum length of trim_ns_left at the
#            5'-end.
#-trim_ns_left <integer>
#            Trim poly-N tail with a minimum length of trim_ns_left at the
#            5'-end.
#-trim_ns_right <integer>
#            Trim poly-N tail with a minimum length of trim_ns_right at the
#            3'-end.
#-trim_qual_left <integer>
#            Trim sequence by quality score from the 5'-end with this
#            threshold score.
#-trim_qual_right <integer>
#            Trim sequence by quality score from the 3'-end with this
#            threshold score.
#-trim_qual_type <string>
#            Type of quality score calculation to use. Allowed options are
#            min, mean, max and sum. [default: min]
#-trim_qual_rule <string>
#            Rule to use to compare quality score to calculated value.
#            Allowed options are lt (less than), gt (greater than) and et
#            (equal to). [default: lt]
#-stats_info
#            Outputs basic information such as number of reads (reads) and
#            total bases (bases).
#-stats_len
#            Outputs minimum (min), maximum (max), range (range), mean
#            (mean), standard deviation (stddev), mode (mode) and mode value
#            (modeval), and median (median) for read length.
#-stats_ns
#            Outputs the number of reads with ambiguous base N (seqswithn),
#            the maximum number of Ns per read (maxn) and the maximum
#            percentage of Ns per read (maxp). The maxn and maxp value are
#            not necessary from the same sequence.
#-stats_all
#            Outputs all available summary statistics.


# Filtering and trimming performed in the following order:
#
# seq_num
# trim_left
# trim_right
# trim_left_p
# trim_right_p
# trim_qual_left
# trim_qual_right
# trim_tail_left
# trim_tail_right
# trim_ns_left
# trim_ns_right
# trim_to_len
# min_len
# max_len
# range_len
# min_qual_score
# max_qual_score
# min_qual_mean
# max_qual_mean
# min_gc
# max_gc
# range_gc
# ns_max_p
# ns_max_n
# noniupac
# lc_method
# derep
# seq_id
# seq_case
# dna_rna
# out_format





#note: runtime (for 2x fq.gz files; 0.5 Gb each): up to 2h30m
#                                               : 14h for largest file



