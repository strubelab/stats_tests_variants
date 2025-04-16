#!/bin/bash -l
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH --job-name=stats
#SBATCH --output=./out_ibex/%x.%j.out
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --mail-user=francisco.guzmanvega@kaust.edu.sa
#SBATCH --mail-type=ALL

conda activate ./env

snakemake --cores 16 --resources mem_mb=30000 --sdm conda

