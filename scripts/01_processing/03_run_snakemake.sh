#!/bin/bash
#SBATCH --job-name=CB_processing                                                              
#SBATCH --time=10:00:00                               
#SBATCH --mem=2G
#SBATCH -n 8 # threaded 
#SBATCH --cpus-per-task=1
#SBATCH -o slurm.CB_processing.job.%j.out
#SBATCH -e slurm.CB_processing.job.%j.err

# activate conda environment
source $HOME/.bash_profile
module load python3
conda activate chemobrain

# change directory to where Snakefile is located
CWD="/tgen_labs/jfryer/kolney/chemobrain/scripts/01_processing"
cd $CWD

# run snakemake
# snakemake -s Snakefile --nolock --latency-wait 15 --executor slurm --profile slurm_profile --jobs 2 --resources mem_mb=16000 time=10:00:00 
# snakemake --nolock -s Snakefile \
#   --jobs 2 \
#   --executor slurm \
#   --profile slurm_profile \
#   --rerun-incomplete \
#   --resources mem_mb=16000 \
#   --default-resources mem_gb=64 walltime=01:30:00 threads=8

snakemake --nolock -s Snakefile --jobs 16 --executor slurm --profile slurm_profile --rerun-incomplete --default-resources mem_mb=64000 ntasks=1 threads=8 runtime=550 cpus_per_task=8
