#!/bin/bash




# Estimate StringTie scalling factors




##### paths and folders


#remember current path
SRCDIR=$(pwd)                                               				

#input file (lists paths to StringTie folders associated to each sample) 
BFL=/home/larshook/LarsH/REVISIT_DC/SCRIPTS/StringTie_files_list.txt

#input file with counts per gene
CPG=/home/larshook/LarsH/REVISIT_DC/NORMALIZED_LIBS/gene_count_matrix.csv

#results folder	
RR=/home/larshook/LarsH/REVISIT_DC/NORMALIZED_LIBS

#output file
OUT_file=StringTie_scalling-factors.txt									

#reference genes
REFG1=leptidea_sinapisG00000003009
REFG2=leptidea_sinapisG00000010515
REFG3=leptidea_sinapisG00000010771



######## read StringTie file list

i=1

while read -r LINE                                                                      
do
   SAMPLE_LIST[$i]=$LINE
   let "i+=1"
done < ${BFL}



### Number of samples
N_SAMPLES=${#SAMPLE_LIST[@]}





#### Compute StringTie scalling factors


## FPKM usually computed using the following formula:   FPKM = (counts * 1e6)/(library_size * gene_length)   [gene_length in kilobases]     
## therefore we have:                                   library_size = (counts * 1e6)/(FPKM * gene_length) 
## and thus:                                            scalling factor = library_size/1e6 = (counts)/(FPKM * gene_length)


echo -e "sample"'\t'"scalling_factor" > $RR/$OUT_file


for i in `seq 1 1 $N_SAMPLES`; do  

   lib_aux=${SAMPLE_LIST[$i]}

   cd $lib_aux/stringtie_2
   
   echo
   echo ${lib_aux: -13}

   FPKM_AUX1="$(grep $REFG1 t_data.ctab | cut -f 12)"
   gene_length_AUX1="$(grep $REFG1 t_data.ctab | cut -f 8)"
   thousand_aux="1000"
   gene_length_AUX1="$(awk '{print $1/$2}' <<<"$gene_length_AUX1 $thousand_aux")"

   FPKM_AUX2="$(grep $REFG2 t_data.ctab | cut -f 12)"
   gene_length_AUX2="$(grep $REFG2 t_data.ctab | cut -f 8)"
   gene_length_AUX2="$(awk '{print $1/$2}' <<<"$gene_length_AUX2 $thousand_aux")"

   FPKM_AUX3="$(grep $REFG3 t_data.ctab | cut -f 12)"
   gene_length_AUX3="$(grep $REFG3 t_data.ctab | cut -f 8)"
   gene_length_AUX3="$(awk '{print $1/$2}' <<<"$gene_length_AUX3 $thousand_aux")"


   Count_aux1a="$(head -1 $CPG | tr -s ',' '\n' | nl -nln |  grep ${lib_aux: -13}| cut -f 1)"        #determine column associated to a specific library
   Count_aux1b="$(grep $REFG1 $CPG | cut -f $Count_aux1a -d ",")"                                   #get counts    
   Count_aux2b="$(grep $REFG2 $CPG | cut -f $Count_aux1a -d ",")"
   Count_aux3b="$(grep $REFG3 $CPG | cut -f $Count_aux1a -d ",")"


   sf_lib_1="$(awk '{print $1/($2*$3)}' <<<"$Count_aux1b $FPKM_AUX1 $gene_length_AUX1")"
   sf_lib_2="$(awk '{print $1/($2*$3)}' <<<"$Count_aux2b $FPKM_AUX2 $gene_length_AUX2")"
   sf_lib_3="$(awk '{print $1/($2*$3)}' <<<"$Count_aux3b $FPKM_AUX3 $gene_length_AUX3")"
   three_aux="3"
   sf_lib_avg="$(awk '{print ($1+$2+$3)/$4}' <<<"$sf_lib_1 $sf_lib_2 $sf_lib_3 $three_aux")"


   #echo $FPKM_AUX1
   #echo $gene_length_AUX1
   #echo $Count_aux1b
   echo $sf_lib_1
   echo $sf_lib_2
   echo $sf_lib_3
   echo $sf_lib_avg

   #sample_aux=${lib_aux#/crex1/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/00_data_mirror_Lars/2e_STAR/}

   echo -e ${lib_aux: -13}'\t'$sf_lib_avg >> $RR/$OUT_file

              
done




