#!/bin/bash -l
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH --job-name=stats
#SBATCH --output=./out_ibex/%x.%j.out
#SBATCH --time=24:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=36
#SBATCH --mail-user=francisco.guzmanvega@kaust.edu.sa
#SBATCH --mail-type=ALL

conda activate ./env

snakemake --cores 36 --resources mem_mb=250000 --sdm conda

