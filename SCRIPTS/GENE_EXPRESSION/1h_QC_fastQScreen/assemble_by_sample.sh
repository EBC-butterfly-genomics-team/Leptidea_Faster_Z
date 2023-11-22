#!/bin/bash -l
#SBATCH -J catting
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -A snic2022-5-34

# load modules
module load bioinfo-tools

#input data files

PATH_MAIN=/home/larshook/LarsH/RNAseq_Lsin_2020       #project folder           
PATH_LOCAL=1x_QC_FINAL			             	             #folder containing input data
N_SAMPLES=28                                                 #number of samples 

# output folder
OUT_FOLDER=1z_QC_concatenated_files



SRCDIR=$(pwd)                                                #remember current path


cd $PATH_MAIN
mkdir $OUT_FOLDER


cd $PATH_MAIN/$PATH_LOCAL

find * -maxdepth 0 -type f >$SRCDIR/list_all_files_c.txt     #gets names of all files present in 
                                                             #input data folder, saves names to file

cd $SRCDIR

grep "\.fastq$" list_all_files_c.txt>input_data_file_names_c.txt   #gets names of all data files (.fastq)


RAW_FN=$SRCDIR/input_data_file_names_c.txt

declare -A READS                                             #declare variable as array

i=1
j=1

while read -r LINE                                 #reads each filename at a time, stores names
                                                   #(paired) in array 
do
   if [ "$j" -eq 1 ]
   then
      READS[$i,1]=$LINE
      let "j=2"
   else
      READS[$i,2]=$LINE
      let "i+=1"
      let "j=1"
   fi

done < ${RAW_FN}



declare -A READS                				#declare variable as array

READS[1,1]=202_S48
READS[1,2]=1                   # lanes 5 and 6

READS[2,1]=203_S49
READS[2,2]=1                   # lanes 5 and 6

READS[3,1]=204_S34
READS[3,2]=0                   # lane 4

READS[4,1]=205_S35
READS[4,2]=0                   # lane 4

READS[5,1]=210_S51
READS[5,2]=1                   # lanes 5 and 6

READS[6,1]=211_S52
READS[6,2]=1                   # lanes 5 and 6

READS[7,1]=212_S36
READS[7,2]=0	               # lane 4

READS[8,1]=213_S53
READS[8,2]=1                   # lanes 5 and 6

READS[9,1]=218_S55
READS[9,2]=1                   # lanes 5 and 6

READS[10,1]=219_S56
READS[10,2]=1                   # lanes 5 and 6

READS[11,1]=220_S57
READS[11,2]=1                   # lanes 5 and 6

READS[12,1]=221_S58
READS[12,2]=1                   # lanes 5 and 6

READS[13,1]=226_S59
READS[13,2]=1                   # lanes 5 and 6

READS[14,1]=227_S38
READS[14,2]=0                   # lane 4

READS[15,1]=233_S60
READS[15,2]=1                   # lanes 5 and 6

READS[16,1]=234_S61
READS[16,2]=1                   # lanes 5 and 6

READS[17,1]=235_S62
READS[17,2]=1                   # lanes 5 and 6



READS[18,1]=201_S47
READS[18,2]=1                   # lanes 5 and 6

READS[19,1]=209_S50
READS[19,2]=1                   # lanes 5 and 6

READS[20,1]=217_S54
READS[20,2]=1                   # lanes 5 and 6


#READS[21,1]=325_S10
#READS[21,2]=0                   # lane 4



# concatenate files associated to the same sample (applies to lane 5 and 6)

for i in `seq 1 1 $N_SAMPLES`; do                

   if [ "${READS[$i,2]}" -eq 1 ]
   then


      cp $PATH_MAIN/$PATH_LOCAL/adjusted-condetri-P5052_${READS[$i,1]}_L005_trim1.tagged_filter.fastq $SNIC_TMP/adjusted-condetri-P5052_${READS[$i,1]}_L005_trim1.tagged_filter.fastq
      cp $PATH_MAIN/$PATH_LOCAL/adjusted-condetri-P5052_${READS[$i,1]}_L005_trim2.tagged_filter.fastq $SNIC_TMP/adjusted-condetri-P5052_${READS[$i,1]}_L005_trim2.tagged_filter.fastq
      cp $PATH_MAIN/$PATH_LOCAL/adjusted-condetri-P5052_${READS[$i,1]}_L006_trim1.tagged_filter.fastq $SNIC_TMP/adjusted-condetri-P5052_${READS[$i,1]}_L006_trim1.tagged_filter.fastq
      cp $PATH_MAIN/$PATH_LOCAL/adjusted-condetri-P5052_${READS[$i,1]}_L006_trim2.tagged_filter.fastq $SNIC_TMP/adjusted-condetri-P5052_${READS[$i,1]}_L006_trim2.tagged_filter.fastq

      cat $SNIC_TMP/adjusted-condetri-P5052_${READS[$i,1]}_L005_trim1.tagged_filter.fastq $SNIC_TMP/adjusted-condetri-P5052_${READS[$i,1]}_L006_trim1.tagged_filter.fastq > $SNIC_TMP/adjusted-condetri-P5052_${READS[$i,1]}_trim1.tagged_filter.fastq
      
      cat $SNIC_TMP/adjusted-condetri-P5052_${READS[$i,1]}_L005_trim2.tagged_filter.fastq $SNIC_TMP/adjusted-condetri-P5052_${READS[$i,1]}_L006_trim2.tagged_filter.fastq > $SNIC_TMP/adjusted-condetri-P5052_${READS[$i,1]}_trim2.tagged_filter.fastq

      gzip -c $SNIC_TMP/adjusted-condetri-P5052_${READS[$i,1]}_trim1.tagged_filter.fastq > $PATH_MAIN/$OUT_FOLDER/adjusted-condetri-P5052_${READS[$i,1]}_trim1.tagged_filter.fq.gz
      
      gzip -c $SNIC_TMP/adjusted-condetri-P5052_${READS[$i,1]}_trim2.tagged_filter.fastq > $PATH_MAIN/$OUT_FOLDER/adjusted-condetri-P5052_${READS[$i,1]}_trim2.tagged_filter.fq.gz

   else

      gzip -c $PATH_MAIN/$PATH_LOCAL/adjusted-condetri-P5052_${READS[$i,1]}_L004_trim1.tagged_filter.fastq > $PATH_MAIN/$OUT_FOLDER/adjusted-condetri-P5052_${READS[$i,1]}_trim1.tagged_filter.fq.gz
      gzip -c $PATH_MAIN/$PATH_LOCAL/adjusted-condetri-P5052_${READS[$i,1]}_L004_trim2.tagged_filter.fastq > $PATH_MAIN/$OUT_FOLDER/adjusted-condetri-P5052_${READS[$i,1]}_trim2.tagged_filter.fq.gz




   fi

done



# runtime: 7h30m


##





