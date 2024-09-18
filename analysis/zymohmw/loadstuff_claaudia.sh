#!/bin/bash
#SBATCH --job-name=basecalling-benchmarks
#SBATCH --output=slurm_logs/%x-%j.out
#SBATCH --error=slurm_logs/%x-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   # Adjust this to the desired number of threads
#SBATCH --mem=16G           # Adjust this to the desired memory allocation
#SBATCH --time=00-10:00:00     # Adjust this to the desired time limit
#SBATCH --mail-user=rhk@bio.aau.dk   # Email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL    # When to send email notifications

module load Mamba/4.14.0-0
source activate /home/bio.aau.dk/ur36rv/.conda/envs/2023-01-basecalling-benchmarks

snakemake --profile profile/