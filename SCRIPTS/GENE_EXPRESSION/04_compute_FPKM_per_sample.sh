#!/bin/bash


## Compute mean FPKM (renormalized) per sex and developmental stage. Redo tables with gene assignment to chromosome (as these files include FPKM values)




################################################

# 1: Path to StringTie main folder containing renormalized libs:

################################################

#Main folder containing StringTie output subfolders for each sample (renormalized)
StringTie=/home/larshook/LarsH/REVISIT_DC/NORMALIZED_LIBS/STAR_norm

#STAR-StringTie outfile name (same for all samples)	
STAR=STAR.Aligned.out.gene_abund.tab		

#OUTPUT FOLDER
RR=/home/larshook/LarsH/REVISIT_DC/NORMALIZED_LIBS

# auxilaiary scripts folder
SRCDIR=/home/larshook/LarsH/REVISIT_DC/SCRIPTS




#######################################

# 2: List with assigned scaffolds

#######################################

ASSIGNED_A_OR_Z=/home/larshook/LarsH/REVISIT_DC/mapped_scaffolds_filtered_above_200_assigned_a_or_z_by_0.95.txt
ASSIGNED_CHR=/home/larshook/LarsH/REVISIT_DC/mapped_scaffolds_filtered_above_200_assigned_to_chromosome_by_0.9.txt





#############################

# 3: Name groups and samples:

#############################

#Groups are combined two-by-two for filtering and comparison (G1 vs G2 etc.)

STAGE1="instar_V"

GROUP1="instar_V_female"
S1=S49
S2=S52
S3=S56

GROUP2="instar_V_male"
S4=S48 
S5=S51
S6=S55

STAGE2="pupa"

GROUP3="pupa_female"
S7=S35
S8=S53
S9=S58

GROUP4="pupa_male"
S10=S34
S11=S36
S12=S57

STAGE3="adult"

GROUP5="adult_female"
S13=S59
S14=S38
S15=S10

GROUP6="adult_male"
S16=S60
S17=S61
S18=S62





############################################

# 4:Run script: ./prepare_stringtie_files.sh

############################################

cd $RR

date > prepare_stringtie_files.log


########################

mkdir -p $RR/samples/$GROUP1

echo
echo preparing $GROUP1 samples...

perl $SRCDIR/1_sort_samples.pl $StringTie/*$S1/stringtie_2/$STAR $RR/samples/$GROUP1/$S1.gene_abund.tab
perl $SRCDIR/1_sort_samples.pl $StringTie/*$S2/stringtie_2/$STAR $RR/samples/$GROUP1/$S2.gene_abund.tab
perl $SRCDIR/1_sort_samples.pl $StringTie/*$S3/stringtie_2/$STAR $RR/samples/$GROUP1/$S3.gene_abund.tab

echo calulating mean FPKM...

cat $RR/samples/$GROUP1/*.tab | perl $SRCDIR/2_calculate_mean.pl $RR/samples/$GROUP1/FPKM_$GROUP1.txt $GROUP1



mkdir -p $RR/samples/$GROUP2

echo preparing $GROUP2 samples...

perl $SRCDIR/1_sort_samples.pl $StringTie/*$S4/stringtie_2/$STAR $RR/samples/$GROUP2/$S4.gene_abund.tab
perl $SRCDIR/1_sort_samples.pl $StringTie/*$S5/stringtie_2/$STAR $RR/samples/$GROUP2/$S5.gene_abund.tab
perl $SRCDIR/1_sort_samples.pl $StringTie/*$S6/stringtie_2/$STAR $RR/samples/$GROUP2/$S6.gene_abund.tab

echo calulating mean FPKM...

cat $RR/samples/$GROUP2/*.tab | perl $SRCDIR/2_calculate_mean.pl $RR/samples/$GROUP2/FPKM_$GROUP2.txt $GROUP2

echo concatenating $STAGE1...
paste $RR/samples/$GROUP1/FPKM_$GROUP1.txt $RR/samples/$GROUP2/FPKM_$GROUP2.txt | cut -f 1-4,8 > $RR/$STAGE1-concatenated.txt

echo assigning genes as A or Z...
perl $SRCDIR/3a_assign_genes.pl $ASSIGNED_A_OR_Z $RR/$STAGE1-concatenated.txt

echo assigning genes to chromosomes...
perl $SRCDIR/3b_assign_genes_to_chromosomes.pl $ASSIGNED_CHR $RR/$STAGE1-concatenated.txt 

echo filtering out non-expressed genes...
perl $SRCDIR/4_filter_genes.pl $RR/$STAGE1-assigned_A_or_Z.txt assigned_A_or_Z $RR/$STAGE1-assigned_A_or_Z-filtered.txt
perl $SRCDIR/4_filter_genes.pl $RR/$STAGE1-assigned_to_chromosomes.txt assigned_to_chromosomes $RR/$STAGE1-assigned_to_chromosomes-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl $RR/$STAGE1-assigned_A_or_Z.txt assigned_A_or_Z female $RR/$STAGE1-assigned_A_or_Z_female-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl $RR/$STAGE1-assigned_A_or_Z.txt assigned_A_or_Z male $RR/$STAGE1-assigned_A_or_Z_male-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl $RR/$STAGE1-assigned_to_chromosomes.txt assigned_to_chromosomes female $RR/$STAGE1-assigned_to_chromosomes_female-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl $RR/$STAGE1-assigned_to_chromosomes.txt assigned_to_chromosomes male $RR/$STAGE1-assigned_to_chromosomes_male-filtered.txt

########################


mkdir -p $RR/samples/$GROUP3

echo
echo preparing $GROUP3 samples...

perl $SRCDIR/1_sort_samples.pl $StringTie/*$S7/stringtie_2/$STAR $RR/samples/$GROUP3/$S7.gene_abund.tab
perl $SRCDIR/1_sort_samples.pl $StringTie/*$S8/stringtie_2/$STAR $RR/samples/$GROUP3/$S8.gene_abund.tab
perl $SRCDIR/1_sort_samples.pl $StringTie/*$S9/stringtie_2/$STAR $RR/samples/$GROUP3/$S9.gene_abund.tab

echo calculating mean FPKM...

cat $RR/samples/$GROUP3/*.tab | perl $SRCDIR/2_calculate_mean.pl $RR/samples/$GROUP3/FPKM_$GROUP3.txt $GROUP3

mkdir -p $RR/samples/$GROUP4

echo preparing $GROUP4 samples...

perl $SRCDIR/1_sort_samples.pl $StringTie/*$S10/stringtie_2/$STAR $RR/samples/$GROUP4/$S10.gene_abund.tab
perl $SRCDIR/1_sort_samples.pl $StringTie/*$S11/stringtie_2/$STAR $RR/samples/$GROUP4/$S11.gene_abund.tab
perl $SRCDIR/1_sort_samples.pl $StringTie/*$S12/stringtie_2/$STAR $RR/samples/$GROUP4/$S12.gene_abund.tab

echo calculating mean FPKM...

cat $RR/samples/$GROUP4/*.tab | perl $SRCDIR/2_calculate_mean.pl $RR/samples/$GROUP4/FPKM_$GROUP4.txt $GROUP4

echo concatenating $STAGE2...
paste $RR/samples/$GROUP3/FPKM_$GROUP3.txt $RR/samples/$GROUP4/FPKM_$GROUP4.txt | cut -f 1-4,8 > $RR/$STAGE2-concatenated.txt

echo assigning genes as A or Z...
perl $SRCDIR/3a_assign_genes.pl $ASSIGNED_A_OR_Z $RR/$STAGE2-concatenated.txt

echo assigning genes to chromosomes...
perl $SRCDIR/3b_assign_genes_to_chromosomes.pl $ASSIGNED_CHR $RR/$STAGE2-concatenated.txt 

echo filtering out non-expressed genes...
perl $SRCDIR/4_filter_genes.pl $RR/$STAGE2-assigned_A_or_Z.txt assigned_A_or_Z $RR/$STAGE2-assigned_A_or_Z-filtered.txt
perl $SRCDIR/4_filter_genes.pl $RR/$STAGE2-assigned_to_chromosomes.txt assigned_to_chromosomes $RR/$STAGE2-assigned_to_chromosomes-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl $RR/$STAGE2-assigned_A_or_Z.txt assigned_A_or_Z female $RR/$STAGE2-assigned_A_or_Z_female-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl $RR/$STAGE2-assigned_A_or_Z.txt assigned_A_or_Z male $RR/$STAGE2-assigned_A_or_Z_male-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl $RR/$STAGE2-assigned_to_chromosomes.txt assigned_to_chromosomes female $RR/$STAGE2-assigned_to_chromosomes_female-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl $RR/$STAGE2-assigned_to_chromosomes.txt assigned_to_chromosomes male $RR/$STAGE2-assigned_to_chromosomes_male-filtered.txt

########################


mkdir -p $RR/samples/$GROUP5

echo
echo preparing $GROUP5 samples...

perl $SRCDIR/1_sort_samples.pl $StringTie/*$S13/stringtie_2/$STAR $RR/samples/$GROUP5/$S13.gene_abund.tab
perl $SRCDIR/1_sort_samples.pl $StringTie/*$S14/stringtie_2/$STAR $RR/samples/$GROUP5/$S14.gene_abund.tab
perl $SRCDIR/1_sort_samples.pl $StringTie/*$S15/stringtie_2/$STAR $RR/samples/$GROUP5/$S15.gene_abund.tab

echo calculating mean FPKM...
cat $RR/samples/$GROUP5/*.tab | perl $SRCDIR/2_calculate_mean.pl $RR/samples/$GROUP5/FPKM_$GROUP5.txt $GROUP5

mkdir -p $RR/samples/$GROUP6

echo preparing $GROUP6 samples...

perl $SRCDIR/1_sort_samples.pl $StringTie/*$S16/stringtie_2/$STAR $RR/samples/$GROUP6/$S16.gene_abund.tab
perl $SRCDIR/1_sort_samples.pl $StringTie/*$S17/stringtie_2/$STAR $RR/samples/$GROUP6/$S17.gene_abund.tab
perl $SRCDIR/1_sort_samples.pl $StringTie/*$S18/stringtie_2/$STAR $RR/samples/$GROUP6/$S18.gene_abund.tab

echo calculating mean FPKM...
cat $RR/samples/$GROUP6/*.tab | perl $SRCDIR/2_calculate_mean.pl $RR/samples/$GROUP6/FPKM_$GROUP6.txt $GROUP6

echo concatenating $STAGE3...
paste $RR/samples/$GROUP5/FPKM_$GROUP5.txt $RR/samples/$GROUP6/FPKM_$GROUP6.txt | cut -f 1-4,8 > $RR/$STAGE3-concatenated.txt

echo assigning genes as A or Z...
perl $SRCDIR/3a_assign_genes.pl $ASSIGNED_A_OR_Z $RR/$STAGE3-concatenated.txt

echo assigning genes to chromosomes...
perl $SRCDIR/3b_assign_genes_to_chromosomes.pl $ASSIGNED_CHR $RR/$STAGE3-concatenated.txt 

echo filtering out non-expressed genes...
perl $SRCDIR/4_filter_genes.pl $RR/$STAGE3-assigned_A_or_Z.txt assigned_A_or_Z $RR/$STAGE3-assigned_A_or_Z-filtered.txt
perl $SRCDIR/4_filter_genes.pl $RR/$STAGE3-assigned_to_chromosomes.txt assigned_to_chromosomes $RR/$STAGE3-assigned_to_chromosomes-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl $RR/$STAGE3-assigned_A_or_Z.txt assigned_A_or_Z female $RR/$STAGE3-assigned_A_or_Z_female-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl $RR/$STAGE3-assigned_A_or_Z.txt assigned_A_or_Z male $RR/$STAGE3-assigned_A_or_Z_male-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl $RR/$STAGE3-assigned_to_chromosomes.txt assigned_to_chromosomes female $RR/$STAGE3-assigned_to_chromosomes_female-filtered.txt
perl $SRCDIR/5_filter_genes_individualSexes.pl $RR/$STAGE3-assigned_to_chromosomes.txt assigned_to_chromosomes male $RR/$STAGE3-assigned_to_chromosomes_male-filtered.txt





#############################################################################
