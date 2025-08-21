#!/bin/bash
#SBATCH --job-name=manhattan_Fcs
#SBATCH --ntasks-per-node=1
#SBATCH -t 1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --error=slurm_logs/manhattan/job.%J.err
#SBATCH --output=slurm_logs/manhattan/job.%J.out
#SBATCH --mem=100G


R CMD BATCH run_plot_manhattan_Fcs_loop.R
