#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 20:00:00
#SBATCH -J samtools_concensus


module load bioinfo-tools

# male
for i in 102_S23 104_S25 106_S27 108_S29

# female
# for i in 103_S24 104_S25 107_S28 108_S29
# do

  in_path=/home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_MAPPING/MARK_DUPLICATES
  in_file=LsinapisSweM-P14458_"$i"-sorted.bam-unique.deduped.bam
  region=LsinapisSweM_CDS_regions.txt
  out_dir=/home/larshook/LarsH/FastZ

  cd "$out_dir"

  rm -f "$i"_concensus_CDS.fa


  while read line;
  do

    /home/larshook/LarsH/SOFTWARE/samtools-1.16.1/samtools consensus -A -r "$line" -o "$out_dir"/"$i"_"$line".fasta "$in_path"/"$in_file"

    # set species id, scaffold and region coords as fasta header and add fasta to multifasta
    awk '/^>/ {gsub(/.fa(sta)?$/,"",FILENAME);printf(">%s\n",FILENAME);next;} {print}' "$i"_"$line".fasta >> "$i"_concensus_CDS.fa

    # remove individual fasta
    rm -f *fasta

  done < "$out_dir"/"$region"

done
