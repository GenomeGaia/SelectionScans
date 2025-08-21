#!/bin/sh
#SBATCH --job-name=q_pbs
#SBATCH --ntasks-per-node=1
#SBATCH -t 0-20:00:00
#SBATCH --cpus-per-task=1
#SBATCH --error=slurm_logs/quantile/job.%J.err
#SBATCH --output=slurm_logs/quantile/job.%J.out
#SBATCH --mem=5G


R CMD BATCH run_quantile_files.R
