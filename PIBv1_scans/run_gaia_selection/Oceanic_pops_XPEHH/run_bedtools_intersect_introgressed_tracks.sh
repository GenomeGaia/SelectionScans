#!/bin/bash
#SBATCH --job-name=bedtools_AI
#SBATCH --ntasks-per-node=1
#SBATCH -t 0-02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --error=slurm_logs/AI/job.%J.err
#SBATCH --output=slurm_logs/AI/job.%J.out
#SBATCH --mem=20G

#SBATCH -a 0-16

#Set the SLURM array task range from 0 to the number of populations to scan minus one

arr=(*0.001.tsv)


i="${SLURM_ARRAY_TASK_ID}"

file_name=${arr[$i]}
name_output1="${file_name//0.001.tsv/0.001_introgressed_corehaps.tsv}"
name_output2="${name_output1//introgressed_corehaps.tsv/introgressed_corehaps_uniquePOP.tsv}"


bedtools intersect -a  <(tail -n+2 ${file_name}) -b <(tail -n+2 ../../key_files/PIBv1_Sprime_putative_adaptive_introgressed_corehaps_OCNnoVanuatu_wGeneList_lumpedpopsrenamed.bed) -wao >  ${name_output1} 

awk 'BEGIN{FS="\t";OFS=FS;}$6==$22{print;}' ${name_output1} | cat ../../key_files/header_file.txt - >  ${name_output2}
