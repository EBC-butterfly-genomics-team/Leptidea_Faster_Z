#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 30:00
#SBATCH -J convert_to_phylip


cd /home/larshook/LarsH/FastZ/DN_DS/FILTERED_ALIGNMENTS

for i in {00000001..00014378}
do

  file=gene_"$i".best.filtered.fas
  
  if [ -f "$file" ];
  then

    # count number of sequences...
    seq_count=$(grep ">" gene_"$i".best.filtered.fas | wc -l)

    # count alignment length...
    ali_len=$(grep -v ">" gene_"$i".best.filtered.fas | wc -c | awk -v seq_count="$seq_count" '{print $1/seq_count-1}')

    # make new file
    printf "$seq_count $ali_len\n" > gene_"$i".best.filtered.phy
    cat gene_"$i".best.filtered.fas | sed 's/>//' >> gene_"$i".best.filtered.phy

  else
  
    echo "$file does not exist."
  
  fi
done
