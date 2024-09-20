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

# Exit on first error and if any variables are unset
set -eu

#source ~/.bashrc  # or ~/.zshrc

# Ensure the PATH is set correctly
#export PATH="/home/bio.aau.dk/ur36rv/.conda/envs/microbench_zymohmw_basecalling-benchmarks/bin:$PATH"

eval "$(conda shell.bash hook)"
conda activate /home/bio.aau.dk/ur36rv/.conda/envs/microbench_zymohmw_basecalling-benchmarks

snakemake --profile profile/

# cluster command not in snakemake 8... downgrading
#snakemake --cluster "sbatch --parsable --output=jobs/{rule}/slurm_%x_%j.out --error=jobs/{rule}/slurm_%x_%j.log --mem={resources.mem_mb} --cpus-per-task={threads} --time 00-10:00:00" -j 400 --use-conda
#snakemake --cluster qsub --jobs 32
#snakemake \
#         --executor cluster-generic \
#         --cluster-generic-submit-cmd 'qsub -N {rule} -q all.q -l h_vmem=8G -pe smp {threads} -V -cwd'