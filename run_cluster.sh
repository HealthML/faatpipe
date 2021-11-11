#!/bin/bash

#SBATCH --job-name=controljob_%j
#SBATCH --output='snakemake.log'
#SBATCH --partition=vcpu
#SBATCH --time=168:00:00
#SBATCH -c 1
#SBATCH --mem 2000

# Initialize conda:
eval "$(conda shell.bash hook)"

snakemake_env="install/snakemake"

if [ ! -d $snakemake_env ]; then
    ./install.sh
fi

conda activate $snakemake_env

snakemake --snakefile Snakefile \
	  "$@" \
          --configfile conf/config.yaml \
	  --profile ./slurm \
          --directory "${PWD}"


