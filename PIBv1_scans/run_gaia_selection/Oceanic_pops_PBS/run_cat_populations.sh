#!/bin/bash

#Hard-coded version:
#populations="Ata Baining-Kagat Baining-Mali Bellona-Rennell Kove Lavongai-Mussau Malaita Mamusi Melamela Nailik-Notsi-Tigak Nakanai-Mangseng Nasioi Santa-Cruz Saposa Sepik-Goroka Tikopia Vella-Lavella"
#for population in $populations; do
#  head -n 1 PIBv1_chr1_Ata_pbs.tsv > ${population}_PIBv1_ALL_chr_pbs.tsv
#  awk 'FNR>1' $(for i in {1..22}; do echo -n "PIBv1_chr${i}_${population}_pbs.tsv "; done) >> ${population}_PIBv1_ALL_chr_pbs.tsv
#done
#Generalized version relying on contents of pops_oceania subdirectory:
while read population; do
  awk 'BEGIN{FS="\t";OFS=FS;filenum=0;}FNR==1{filenum+=1;}filenum==1&&FNR==1{print;}FNR>1{print;}' PIBv1_chr{1..22}_${population}_pbs.tsv > ${population}_PIBv1_ALL_chr_pbs.tsv
done < <(find pops_oceania/ -name "pop1_*.txt" -print | sed 's/pops_oceania\/pop1_//;s/.txt//;' | sort)
