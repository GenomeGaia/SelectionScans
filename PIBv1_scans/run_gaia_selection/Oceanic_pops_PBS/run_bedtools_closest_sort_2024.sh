#!/bin/bash
#SBATCH --job-name=bedtools
#SBATCH --ntasks-per-node=1
#SBATCH -t 1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --error=slurm_logs/bclosest/job.%J.err
#SBATCH --output=slurm_logs/bclosest/job.%J.out
#SBATCH --mem=20G
#SBATCH -a 0-16

#Set the SLURM array task range from 0 to the number of populations to scan minus one

arr=(*.tsv)

i="${SLURM_ARRAY_TASK_ID}"

file_name="${arr[$i]}"

output1="${file_name//.tsv/_sort.bed}"
output2="${output1/_sort.bed/_sort_closest.tsv}"

#
sort -k1,1V -k2,2n -k3,3n "${file_name}" | sed '$ d' > "${output1}"
#
bedtools closest -a "${output1}" -b "../../key_files/gencode.v38lift37.annotation.codingGENES_sorted.bed" -d -t all -g "../../key_files/hs37d5_wChr.genome" > "$(basename ${output2})"

