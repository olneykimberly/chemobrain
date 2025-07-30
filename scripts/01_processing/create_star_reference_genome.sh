#!/bin/bash
#SBATCH --job-name=STAR_mouse_ref                                                              
#SBATCH --time=05:00:00                               
#SBATCH --mem=64G
#SBATCH -n 1 # threaded 
#SBATCH --cpus-per-task=8
#SBATCH -o slurm.STAR_mouse_ref.out
#SBATCH -e slurm.STAR_mouse_ref.err

# activate conda environment
source $HOME/.bash_profile
module load python3
conda activate chemobrain

# change directory to where Snakefile is located
CWD="/tgen_labs/jfryer/projects/references/mouse"
cd $CWD

STAR --runMode genomeGenerate --runThreadN 8 --genomeDir refdata-gex-GRCm39-2024-A_STARv2.7.11_150sjdb  --genomeFastaFiles /tgen_labs/jfryer/projects/references/mouse/refdata-gex-GRCm39-2024-A/fasta/genome.fa --sjdbGTFfile /tgen_labs/jfryer/projects/references/mouse/refdata-gex-GRCm39-2024-A/genes/genes.gtf --sjdbOverhang 150
