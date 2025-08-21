#!/bin/sh
#SBATCH --job-name=final_pbs
#SBATCH --ntasks-per-node=1
#SBATCH -t 1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --error=slurm_logs/final_pbs/job.%J.err
#SBATCH --output=slurm_logs/final_pbs/job.%J.out
#SBATCH --mem=20G

#SBATCH -a 0-16

#Set the SLURM array task range from 0 to the number of populations to scan minus one

./run_topGENES_quantiles_final_2024.py
