#!/bin/sh
#SBATCH --job-name=final-stats
#SBATCH --ntasks-per-node=1
#SBATCH -t 1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --error=slurm_logs/windows_pbs-xpehh/job.%J.err
#SBATCH --output=slurm_logs/windows_pbs-xpehh/job.%J.out
#SBATCH --mem=50G

#SBATCH -a 0-16

#Set the SLURM array task range from 0 to the number of populations to scan minus one

/usr/bin/time -v ./run_final-stats_HMP_pbs-xpehh_2024.py
