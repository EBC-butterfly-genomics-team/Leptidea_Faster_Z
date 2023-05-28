#!/bin/bash -l

#SBATCH -A naiss2023-5-52
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -J calculate_mean_FPKM


rna_path=/home/larshook/LarsH/FastZ/RNAseq_Lsin/2e_STAR
out_path=/home/larshook/LarsH/FastZ/RNAseq_Lsin



rm -f $out_path/unexpressed_genes_list.txt 
rm -f $out_path/*mean_FPKM_expressed_genes.txt
rm -f $out_path/*mean_FPKM_female_expressed_genes.txt
rm -f $out_path/*mean_FPKM_male_expressed_genes.txt
rm -f $out_path/classified_expressed_genes.txt


# calculate means...

for i in {00000001..00014378}
do

    # larva female
    larva_female_1=$(grep -w "LSSWEG$i" $rna_path/P5052_203_S49/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')
    larva_female_2=$(grep -w "LSSWEG$i" $rna_path/P5052_211_S52/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')
    larva_female_3=$(grep -w "LSSWEG$i" $rna_path/P5052_219_S56/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')

    mean_larva_female=$(echo $larva_female_1 $larva_female_2 $larva_female_3 | awk '{print ($1+$2+$3)/3}')

    # larva male
    larva_male_1=$(grep -w "LSSWEG$i" $rna_path/P5052_202_S48/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')
    larva_male_2=$(grep -w "LSSWEG$i" $rna_path/P5052_210_S51/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')
    larva_male_3=$(grep -w "LSSWEG$i" $rna_path/P5052_218_S55/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')

    mean_larva_male=$(echo $larva_male_1 $larva_male_2 $larva_male_3 | awk '{print ($1+$2+$3)/3}')

    mean_larva_total=$(echo $larva_female_1 $larva_female_2 $larva_female_3 $larva_male_1 $larva_male_2 $larva_male_3 | awk '{print ($1+$2+$3+$4+$5+$6)/6}')


    # pupa female
    pupa_female_1=$(grep -w "LSSWEG$i" $rna_path/P5052_205_S35/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')
    pupa_female_2=$(grep -w "LSSWEG$i" $rna_path/P5052_213_S53/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')
    pupa_female_3=$(grep -w "LSSWEG$i" $rna_path/P5052_221_S58/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')

    mean_pupa_female=$(echo $pupa_female_1 $pupa_female_2 $pupa_female_3 | awk '{print ($1+$2+$3)/3}')

    # pupa male
    pupa_male_1=$(grep -w "LSSWEG$i" $rna_path/P5052_204_S34/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')
    pupa_male_2=$(grep -w "LSSWEG$i" $rna_path/P5052_212_S36/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')
    pupa_male_3=$(grep -w "LSSWEG$i" $rna_path/P5052_220_S57/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')

    mean_pupa_male=$(echo $pupa_male_1 $pupa_male_2 $pupa_male_3 | awk '{print ($1+$2+$3)/3}')

    mean_pupa_total=$(echo $pupa_female_1 $pupa_female_2 $pupa_female_3 $pupa_male_1 $pupa_male_2 $pupa_male_3 | awk '{print ($1+$2+$3+$4+$5+$6)/6}')


    # adult female
    adult_female_1=$(grep -w "LSSWEG$i" $rna_path/P7553_325_S10/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')
    adult_female_2=$(grep -w "LSSWEG$i" $rna_path/P5052_227_S38/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')
    adult_female_3=$(grep -w "LSSWEG$i" $rna_path/P5052_226_S59/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')

    mean_adult_female=$(echo $adult_female_1 $adult_female_2 $adult_female_3 | awk '{print ($1+$2+$3)/3}')

    # adult male
    adult_male_1=$(grep -w "LSSWEG$i" $rna_path/P5052_233_S60/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')
    adult_male_2=$(grep -w "LSSWEG$i" $rna_path/P5052_234_S61/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')
    adult_male_3=$(grep -w "LSSWEG$i" $rna_path/P5052_235_S62/stringtie_2/STAR.Aligned.out.gene_abund.tab | awk '{print $8}')

    mean_adult_male=$(echo $adult_male_1 $adult_male_2 $adult_male_3 | awk '{print ($1+$2+$3)/3}')
    
    mean_adult_total=$(echo $adult_female_1 $adult_female_2 $adult_female_3 $adult_male_1 $adult_male_2 $adult_male_3 | awk '{print ($1+$2+$3+$4+$5+$6)/6}')




    # print unexpressed genes to one file and expressed to others...

  
      if [[ $mean_larva_total = "0" ]]
      then
        printf "LSSWEG$i\tnon\ttotal\tlarva\t"$mean_larva_total"\n" >> $out_path/unexpressed_genes_list.txt
      else
        printf "LSSWEG$i\t"$mean_larva_total"\n" >> $out_path/larva_mean_FPKM_expressed_genes.txt
      fi

      if [[ "$mean_larva_female" = "0" ]]
      then
        printf "LSSWEG$i\tnon\tfemale\tlarva\t"$mean_larva_female"\n" >> $out_path/unexpressed_genes_list.txt
      else
        printf "LSSWEG$i\t"$mean_larva_female"\n" >> $out_path/larva_mean_FPKM_female_expressed_genes.txt
      fi

      if [[ "$mean_larva_male" = "0" ]]
      then
        printf "LSSWEG$i\tnon\tmale\tlarva\t"$mean_larva_male"\n" >> $out_path/unexpressed_genes_list.txt
      else
        printf "LSSWEG$i\t"$mean_larva_male"\n" >> $out_path/larva_mean_FPKM_male_expressed_genes.txt
      fi


      if [[ $mean_pupa_total = "0" ]]
      then
        printf "LSSWEG$i\tnon\ttotal\tpupa\t"$mean_pupa_total"\n" >> $out_path/unexpressed_genes_list.txt
      else
	printf "LSSWEG$i\t"$mean_pupa_total"\n" >> $out_path/pupa_mean_FPKM_expressed_genes.txt
      fi

      if [[ "$mean_pupa_female" = "0" ]]
      then
	printf "LSSWEG$i\tnon\tfemale\tpupa\t"$mean_pupa_female"\n" >> $out_path/unexpressed_genes_list.txt
      else
	printf "LSSWEG$i\t"$mean_pupa_female"\n" >> $out_path/pupa_mean_FPKM_female_expressed_genes.txt
      fi

      if [[ "$mean_pupa_male" = "0" ]]
      then
	printf "LSSWEG$i\tnon\tmale\tpupa\t"$mean_pupa_male"\n" >> $out_path/unexpressed_genes_list.txt
      else
	printf "LSSWEG$i\t"$mean_pupa_male"\n" >> $out_path/pupa_mean_FPKM_male_expressed_genes.txt
      fi


      if [[ $mean_adult_total = "0" ]]
      then
	printf "LSSWEG$i\tnon\ttotal\tadult\t"$mean_adult_total"\n" >> $out_path/unexpressed_genes_list.txt
      else
	printf "LSSWEG$i\t"$mean_adult_total"\n" >> $out_path/adult_mean_FPKM_expressed_genes.txt
      fi

      if [[ "$mean_adult_female" = "0" ]]
      then
        printf "LSSWEG$i\tnon\tfemale\tadult\t"$mean_adult_female"\n" >> $out_path/unexpressed_genes_list.txt
      else
        printf "LSSWEG$i\t"$mean_adult_female"\n" >> $out_path/adult_mean_FPKM_female_expressed_genes.txt
      fi

      if [[ "$mean_adult_male" = "0" ]]
      then
        printf "LSSWEG$i\tnon\tmale\tadult\t"$mean_adult_male"\n" >> $out_path/unexpressed_genes_list.txt
      else
	printf "LSSWEG$i\t"$mean_adult_male"\n" >> $out_path/adult_mean_FPKM_male_expressed_genes.txt
      fi



done



# sort and classify genes as low, mid, high...

for stage in larva pupa adult
do
  total_group_size=$(wc -l $out_path/"$stage"_mean_FPKM_expressed_genes.txt | awk '{printf("%.0f", $1/3)}')
  
  sort -V -k2,2 $out_path/"$stage"_mean_FPKM_expressed_genes.txt | awk -v size=$total_group_size -v stage=$stage -v OFS="\t" 'NR>=1 && NR<=size {print $1, "low", "total", stage}; \
	NR>size && NR<=2*size {print $1, "mid", "total", stage}; \
	NR>2*size {print $1, "high", "total", stage}' >> $out_path/classified_expressed_genes.txt


  female_group_size=$(wc -l $out_path/"$stage"_mean_FPKM_female_expressed_genes.txt | awk '{printf("%.0f", $1/3)}')

  sort -V -k2,2 $out_path/"$stage"_mean_FPKM_female_expressed_genes.txt | awk -v size=$female_group_size -v stage=$stage -v OFS="\t" 'NR>=1 && NR<=size {print $1, "low", "female", stage}; \
        NR>size && NR<=2*size {print $1, "mid", "female", stage}; \
        NR>2*size {print $1, "high", "female", stage}' >> $out_path/classified_expressed_genes.txt


  male_group_size=$(wc -l $out_path/"$stage"_mean_FPKM_male_expressed_genes.txt | awk '{printf("%.0f", $1/3)}')

  sort -V -k2,2 $out_path/"$stage"_mean_FPKM_male_expressed_genes.txt | awk -v size=$female_group_size -v stage=$stage -v OFS="\t" 'NR>=1 && NR<=size {print $1, "low", "male", stage}; \
        NR>size && NR<=2*size {print $1, "mid", "male", stage}; \
        NR>2*size {print $1, "high", "male", stage}' >> $out_path/classified_expressed_genes.txt

done

# combine output...

cat $out_path/unexpressed_genes_list.txt >> $out_path/classified_expressed_genes.txt
