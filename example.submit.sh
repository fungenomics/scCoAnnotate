#!/bin/sh
#SBATCH --job-name=scCoAnnotate
#SBATCH --account=rrg-kleinman 
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=60GB 

module load scCoAnnotate/2.0

# path to snakefile and config 
snakefile=<path to snakefile>
config=<path to configfile>

# unlock directory incase of previous errors
snakemake -s ${snakefile} --configfile ${config} --unlock 

# run workflow 
snakemake -s ${snakefile} --configfile ${config} --cores 5

