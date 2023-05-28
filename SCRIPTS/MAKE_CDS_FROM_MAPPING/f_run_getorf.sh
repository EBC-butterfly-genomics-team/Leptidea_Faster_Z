#!/bin/bash -l


module load bioinfo-tools
module load emboss
module load bbmap


# trim incomplete annotations that don't have start/stop codons...


main_path=/home/larshook/LarsH/FastZ/DN_DS

cd "$main_path"
mkdir -p ORFS


for species in LSSWET 102_S23_ 104_S25_ 106_S27_ 108_S29_
do

  for i in {00000001..00014378}
  do

    getorf \
	-sequence CDS/"$species""$i".fasta \
        -find 3 \
	-reverse N \
        -outseq "$species""$i"_temp.fasta

    # select longest orf and reformat header...

    sortbyname.sh in="$species""$i"_temp.fasta out=stdout.fasta length descending |\
    	awk '/^>/{if(N)exit;++N;} {print;}' |\
    		sed 's/\(.*\)_/\1 /' |\
			awk '{print $1}' > "$main_path"/ORFS/"$species""$i".fasta

  rm -f "$species""$i"_temp.fasta

  done
done
