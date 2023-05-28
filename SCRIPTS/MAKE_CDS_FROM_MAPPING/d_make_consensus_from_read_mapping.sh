#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J samtools_concensus


# make cds consensus sequences for the outgroups based on read mapping against L. sinapis...

module load bioinfo-tools


for i in 102_S23 104_S25 106_S27 108_S29
do

  in_path=/home/larshook/LarsH/FastZ/LEPTIDEA_OUTGROUPS/READ_MAPPING/MARK_DUPLICATES
  in_file=LsinapisSweM-P14458_"$i"-sorted.bam-unique.deduped.bam
  region=LsinapisSweM_CDS_regions.txt
  out_dir=/home/larshook/LarsH/FastZ/DN_DS

  cd "$out_dir"
  rm -f "$i"_CDS.fa


  while read line;
  do

    /home/larshook/LarsH/SOFTWARE/samtools-1.16.1/samtools consensus -A -r "$line" -o "$out_dir"/"$i"_"$line".fasta "$in_path"/"$in_file"

    # set species id, scaffold and region coords as fasta header and add fasta to multifasta
    awk '/^>/ {gsub(/.fa(sta)?$/,"",FILENAME);printf(">%s\n",FILENAME);next;} {print}' "$i"_"$line".fasta >> "$i"_CDS.fa

    # remove individual fasta
    rm -f "$i"_"$line".fasta

  done < "$out_dir"/"$region"
done
