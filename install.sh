#!/bin/bash

snakemake_env="install/snakemake"
set -e

# Install the snakemake environment
if [ ! -d $snakemake_env ]; then
    >&2 echo "Installing snakemake conda environment in $snakemake_env ..." 
    mkdir -p install
    conda env create -f env/snakemake.yml --prefix $snakemake_env
fi

if [ ! -d ./bin ]; then
    mkdir bin
fi

# Unzip some data that is shipped with the pipeline
if [ -f data/regions/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.rm_chr.bed.gz ]; then
    gunzip data/regions/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.rm_chr.bed.gz
fi

if [ -f data/regions/xgen_plus_spikein.b38.bed.gz ]; then
    gunzip data/regions/xgen_plus_spikein.b38.bed.gz
fi

# Install Plink 2.0
if [ ! -f ./bin/plink2 ]; then
    >&2 echo "Downloading Plink 2.0 binaries"
    cd bin
    wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_avx2.zip
    unzip plink2_linux_avx2.zip && rm plink2_linux_avx2.zip
    cd ..
fi

# Install Plink 1.9
if [ ! -f ./bin/plink ]; then
    >&2 echo "Downloading Plink 1.9 binaries"
    cd bin
    wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip
    unzip plink_linux_x86_64_20201019.zip && rm plink_linux_x86_64_20201019.zip
    cd ..
fi
    

