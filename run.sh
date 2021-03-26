#!/bin/bash

snakemake_env="install/snakemake"

if [ ! -d $snakemake_env ]; then
    ./install.sh
fi

source activate $snakemake_env

snakemake --snakefile Snakefile \
          --use-conda \
          --keep-going \
          --configfile conf/config.yaml \
          --directory "${PWD}" \
          --cores 36 \
          "${@}"



