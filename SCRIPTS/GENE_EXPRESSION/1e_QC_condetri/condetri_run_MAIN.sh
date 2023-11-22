###### START UPPMAX RUNS (CALLS condetri_run_query.sh BASH FILE FOR EACH SET OF PAIRED-END READS)

echo
echo "Starting Uppmax jobs ..."
echo



#input data files

#PATH_MAIN=/proj/uppstore2017185/b2014034/nobackup/Lars/RNAseq_Lsin_2020    #project folder           
PATH_MAIN=/home/larshook/LarsH/RNAseq_Lsin_2020
PATH_LOCAL=1d_QC_filterPolyA_prinseq                         #folder containing input data
N_SAMPLES=42                                                 #number of samples 

SRCDIR=$(pwd)                                                #remember current path

cd $PATH_MAIN/$PATH_LOCAL

rm -f *singletons.fastq     #CAUTION!! -> deletes all files with singleton reads (not needed)
rm -f prinseq_REJECTED*     #          -> deletes all files with rejected reads

find * -maxdepth 0 -type f >$SRCDIR/list_all_files.txt       #gets names of all files present in 
                                                             #input data folder, saves names to file

cd $SRCDIR

grep "\.fastq$" list_all_files.txt>input_data_file_names.txt   #gets names of all data files (.fastq)


RAW_FN=$SRCDIR/input_data_file_names.txt

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
   echo ${READS[$i,$j]}
done < ${RAW_FN}




# prepare output folder

SRCDIR=$(pwd)                                    #remember path

cd $PATH_MAIN/                                   #cd to project folder
OUT_FOLDER=1e_QC_condetri/                       #main output folder

if [ -d "$OUT_FOLDER" ]                          #check whether main output folder already exists
   then
      { echo "Error: $PATH_MAIN/$OUT_FOLDER folder already exists."; echo ; exit 1; }
   else
      mkdir $OUT_FOLDER                          #create output folder
fi



# start jobs (60 runs, one for each sample/lane)

cd $SRCDIR

for i in `seq 1 1 $N_SAMPLES`; do                #loop starting all jobs
   RNAFILE_A=${READS[$i,1]}
   RNAFILE_B=${READS[$i,2]}
   echo $RNAFILE_A                               #display fastq file name
   echo $RNAFILE_B
   sbatch condetri_run_query.sh $PATH_MAIN/$PATH_LOCAL/$RNAFILE_A \
                                $PATH_MAIN/$PATH_LOCAL/$RNAFILE_B \
                                $RNAFILE_A \
                                $RNAFILE_B \
                                $PATH_MAIN/$OUT_FOLDER       #starts job
   sleep 1                                       #pauses for 1 sec
   echo
done

echo
jobinfo -u $USER                                    #check job status
echo
