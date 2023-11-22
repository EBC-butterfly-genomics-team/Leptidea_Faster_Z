#!/bin/bash -l
#SBATCH -J adjust_pf
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 01:00:00
#SBATCH -A snic2022-5-34


# load modules
module load bioinfo-tools
module load python/3.5.0


#input data files (DM)

PATH_MAIN=/home/larshook/LarsH/RNAseq_Lsin_2020       #project folder
PATH_LOCAL=1h_QC_fastQScreen                                 #folder containing input data
N_SAMPLES=1                                                 #number of samples

SRCDIR=$(pwd)                                                #remember current path

cd $PATH_MAIN/$PATH_LOCAL

#find * -maxdepth 0 -type f >$SRCDIR/list_all_files2.txt      #gets names of all files present in 
                                                             #input data folder, saves names to file

cd $SRCDIR

#grep "\filter.fastq$" list_all_files2.txt>input_data_file_names2.txt   #gets names of all data files (filter.fastq)


RAW_FN=$SRCDIR/input_data_file_names2.txt



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
#   echo ${READS[$i,$j]}
done < ${RAW_FN}





# run jobs

#SRCDIR=$(pwd)                                      #remember current path
cd $PATH_MAIN/$PATH_LOCAL/                         #cd to output folder
cp $SRCDIR/adjust_paired_files.py .                #copy python script to output folder

for i in `seq 1 1 $N_SAMPLES`; do                  #loop starting all jobs
   python3 adjust_paired_files.py ${READS[$i,1]} ${READS[$i,2]}
done

rm $PATH_MAIN/$PATH_LOCAL/adjust_paired_files.py   #delete copy of python script

#runtime:   up to 10 min per paired-end (4h15 for all 120 files)  
#note:		need at least 4 cores as memory requirements are around 30GB


