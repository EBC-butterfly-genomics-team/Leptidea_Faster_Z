#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J remove_truncated_alignments

# remove alignments if > 50% of initial alignment is trimmed...

main_path=/home/larshook/LarsH/FastZ/DN_DS

cd "$main_path"

mkdir -p TRUNCATED_ALIGNMENTS


for i in {00000001..00014378}
do

  # get sizes pre and post alignment...

  pre_size=$(wc -c "$main_path"/ALIGNMENTS/gene_"$i".best.fas | awk '{print $1}')
  echo "$pre_size"

  post_size=$(wc -c "$main_path"/FILTERED_ALIGNMENTS/gene_"$i".best.filtered.fas | awk '{print $1}')
  echo "$post_size"

  pp_frac=$(echo "scale=2; "$post_size" / "$pre_size"" | bc )

  # remove alignment if < 50%

  echo "$pp_frac"

  remove_align=$(echo "$pp_frac" | awk '{if ($1<0.50) print "yes"}')

  if [[ $remove_align = "yes" ]] 
  then
    mv "$main_path"/FILTERED_ALIGNMENTS/gene_"$i".best.filtered.fas "$main_path"/TRUNCATED_ALIGNMENTS/
  fi

done
