###### START UPPMAX RUNS 

# Performs read mapping with STAR and transcript quantification with StringTie and Cufflinks

echo
echo "Starting Uppmax jobs ..."
echo

#input data files

PATH_MAIN=/home/larshook/LarsH/FastZ/RNAseq_Lsin		#project folder

									#folder containing input data
PATH_LOCAL=1z_QC_concatenated_files

N_SAMPLES=21                                        			#number of samples


SRCDIR=$(pwd)                                               #remember current path



### get file names


cd $PATH_MAIN/$PATH_LOCAL

find * -maxdepth 0 -type f >$SRCDIR/list_all_files.txt       #gets names of all files present in 
                                                             #input data folder, saves names to file

cd $SRCDIR

grep "\.fastq.gz$" list_all_files.txt>input_data_file_names.txt   #gets names of all data files (.fastq)
grep "\.fq.gz$" list_all_files.txt>input_data_file_names.txt

RAW_FN=$SRCDIR/input_data_file_names.txt

declare -A READS                                   #declare variable as array

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








# create main output folder

cd $PATH_MAIN/                                     #cd to project folder
#OUT_FOLDER=2e_STAR				                #main output folder
#OUT_FOLDER=2e_STAR_old
OUT_FOLDER=2e_STAR_FEMALE

if [ -d "$OUT_FOLDER" ]                            #check whether main output folder already exists
   then
      { echo "Error: $PATH_MAIN/$OUT_FOLDER folder already exists."; echo ; exit 1; }
   else
      mkdir $OUT_FOLDER                            #create main output folder
fi



#start jobs

for i in `seq 1 1 $N_SAMPLES`; do                  

   cd $PATH_MAIN/$OUT_FOLDER                    	#cd to project folder

   RNAFILE_A=${READS[$i,1]}
   RNAFILE_B=${READS[$i,2]}

   AUX_DIR=${READS[$i,1]}
   RNADIR=${AUX_DIR#adjusted-condetri-}
   OUT_FOLDER2=${RNADIR%_trim1.tagged_filter.fq.gz}       #output folder
   mkdir ${OUT_FOLDER2}                         	#create output folder for each sample
   cd $PATH_MAIN/$OUT_FOLDER/${OUT_FOLDER2}

   cd $SRCDIR

   echo $RNAFILE_A                              #display fastq file names (paired)
   echo $RNAFILE_B

   sbatch $SRCDIR/STAR_run_query.sh \
                           $PATH_MAIN/$PATH_LOCAL/$RNAFILE_A \
                           $PATH_MAIN/$PATH_LOCAL/$RNAFILE_B \
                           $RNAFILE_A \
                           $RNAFILE_B \
                           $PATH_MAIN/$OUT_FOLDER/${OUT_FOLDER2} 
   
   sleep 1                                                  #pauses for 1 sec
   echo

done


squeue -u $USER                                    #check job status
echo
