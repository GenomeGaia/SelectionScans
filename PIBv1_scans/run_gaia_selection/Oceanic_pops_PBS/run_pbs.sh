#!/bin/sh
#SBATCH --job-name=pbs
#SBATCH --ntasks-per-node=1
#SBATCH -t 1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --error=slurm_logs/pbs/job_%j.stderr
#SBATCH --output=slurm_logs/pbs/job_%j.stdout
#SBATCH --mem=20G

#SBATCH -a 0-373

#Set the SLURM array task range from 0 to the the number of autosomes times the number of populations to scan minus one
#e.g. 22 autosomes * 17 populations to scan = 374, 374 - 1 = 373, therefore -a 0-373

/usr/bin/time -v ./run_pbs.py
