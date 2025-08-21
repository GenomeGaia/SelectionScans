#!/bin/bash
#SBATCH --job-name=fcs
#SBATCH --ntasks-per-node=1
#SBATCH -t 1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --error=slurm_logs/fcs/job.%J.err
#SBATCH --output=slurm_logs/fcs/job.%J.out
#SBATCH --mem=50G


R CMD BATCH run_FCS.R
