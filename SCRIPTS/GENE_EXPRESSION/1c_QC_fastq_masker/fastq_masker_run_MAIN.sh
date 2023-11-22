###### START UPPMAX RUNS 

echo
echo "Starting Uppmax jobs ..."
echo

#input data files

PATH_MAIN=/proj/uppstore2017185/b2014034/nobackup/Lars/RNAseq_Lsin_2020           #project folder
PATH_LOCAL=1b_QC_TrimGalore                                      #folder containing input data
N_SAMPLES=1                                                      #number of samples

SRCDIR=$(pwd)                                                    #remember current path

cd $PATH_MAIN/$PATH_LOCAL

#find * -maxdepth 0 -type f >$SRCDIR/list_all_files.txt           #gets names of all files present in input data folder, saves names to file

cd $SRCDIR

#grep ".fq.gz" list_all_files.txt>input_data_file_names.txt       #names of all data files (fq.gz)

RAW_FN=$SRCDIR/input_data_file_names.txt

declare -A READS                                   #declare variable as array

i=1
j=1

while read -r LINE                                 #reads each filename at a time, stores names (paired) in array 
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
   echo ${READS[$i,$j]}
done < ${RAW_FN}




# prepare output folder

SRCDIR=$(pwd)                                      #remember current path

cd $PATH_MAIN/                                     #cd to project folder
OUT_FOLDER=1c_QC_fastq_masker/                     #output folder

#if [ -d "$OUT_FOLDER" ]                            #check whether output folder already exists
#   then
#      { echo "Error: $PATH_MAIN/$OUT_FOLDER folder already exists."; echo ; exit 1; }
#   else
#      mkdir $OUT_FOLDER                            #create output folder
#fi



# start jobs (60 runs, one for each sample/lane)

cd $SRCDIR

for i in `seq 1 1 $N_SAMPLES`; do                  #loop starting all jobs
   RNAFILE_A=${READS[$i,1]}
   RNAFILE_B=${READS[$i,2]}
   echo $RNAFILE_A                                 #display fastq file names (paired)
   echo $RNAFILE_B
   sbatch fastq_masker_run_query.sh $PATH_MAIN/$PATH_LOCAL/$RNAFILE_A \
                                    $PATH_MAIN/$PATH_LOCAL/$RNAFILE_B \
                                    $RNAFILE_A \
                                    $RNAFILE_B \
                                    $PATH_MAIN/$OUT_FOLDER/                                    
   sleep 1                                                  #pauses for 1 sec
   echo
done


jobinfo -u $USER                                    #check job status
echo
