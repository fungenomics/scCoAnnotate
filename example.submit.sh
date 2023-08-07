#!/bin/sh
#SBATCH --job-name=scCoAnnotate
#SBATCH --account=rrg-kleinman 
#SBATCH --output=logs/%x.out
#SBATCH --error=logs/%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=
#SBATCH --time=
#SBATCH --mem-per-cpu=

module load scCoAnnotate/2.0

# path to snakefile and config 
snakefile=""
config=""

# unlock directory incase of previous errors
snakemake -s ${snakefile} --configfile ${config} --unlock 

# run workflow 
snakemake -s ${snakefile} --configfile ${config} --cores 2

