#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J get_orf

module load bioinfo-tools
module load emboss
module load bbmap


main_path=/home/larshook/LarsH/FastZ/DN_DS

cd "$main_path"
mkdir -p ORFS


for species in 102_S23 104_S25 106_S27 108_S29 LSSWET
do

  for i in {00000001..00014378}
  do

    getorf \
	-sequence "$species"_"$i".fasta \
        -find 3 \
	-reverse N \
        -outseq "$species"_"$i"_temp.fasta

    # select longest orf and reformat header...

    sortbyname.sh in="$species"_"$i"_temp.fasta out=stdout.fasta length descending |\
    	awk '/^>/{if(N)exit;++N;} {print;}' |\
    		sed 's/\(.*\)_/\1 /' |\
			awk '{print $1}' > "$main_path"/ORFS/"$species"_"$i".fasta

  rm -f "$species"_"$i"_temp.fasta

  done
done
