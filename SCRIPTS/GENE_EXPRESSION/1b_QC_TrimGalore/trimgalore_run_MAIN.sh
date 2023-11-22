###### START UPPMAX RUNS (CALLS trimgalore_run_query.sh BASH FILE FOR EACH SET OF PAIRED-END READS)

echo
echo "Starting Uppmax jobs ..."
echo

#input data files

PATH_RAW=/proj/uppstore2017185/b2014034/private/raw_data/RNAseq        						     #path to where raw files are located (top folder)
PATH_MAIN=/proj/uppstore2017185/b2014034/nobackup/Lars/RNAseq_Lsin_2020             					  #path to project folder (my analysis)
RAW_FN=/proj/uppstore2017185/b2014034/nobackup/Lars/RNAseq_Lsin_2020/SCRIPTS/1b_QC_TrimGalore/raw_data_file_names.txt        #file with names of files containing raw reads

N_SAMPLES=42                                       #number of samples

declare -A READS                                   #declare variable as array

i=1
j=1

while read -r LINE                                 #reads each filename at a time, stores names (paired) in array 
do
   if [ "$j" -eq 1 ]
   then
      READS[$i,1]=$PATH_RAW/$LINE
      let "j=2"
   else
      READS[$i,2]=$PATH_RAW/$LINE
      let "i+=1"
      let "j=1"
   fi
   echo ${READS[$i,$j]}
done < ${RAW_FN}



# prepare output folder

SRCDIR=$(pwd)                                      #remember current path

cd $PATH_MAIN/                                     #cd to project folder
OUT_FOLDER=1b_QC_TrimGalore/                       #output folder

if [ -d "$OUT_FOLDER" ]                            #check whether output folder already exists
   then
      { echo "Error: $PATH_MAIN/$OUT_FOLDER folder already exists."; echo ; exit 1; }
   else
      mkdir $OUT_FOLDER                            #create output folder
fi



# start jobs

cd $SRCDIR

for i in `seq 1 1 $N_SAMPLES`; do                  
   RNAFILE_A=${READS[$i,1]}
   RNAFILE_B=${READS[$i,2]}
   echo $RNAFILE_A                                          #display fastq file names (paired)
   echo $RNAFILE_B
   sbatch trimgalore_run_query.sh $RNAFILE_A \
                                  $RNAFILE_B \
                                  $PATH_MAIN/$OUT_FOLDER/   #starts job
   sleep 1                                                  #pauses for 1 sec
   echo
done


jobinfo -u $USER                                    #check job status
echo
