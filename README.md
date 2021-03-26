# Functional Annotation and Association Testing Pipeline
Pipeline to run variant effect prediction and rare-variant association tests for continuous phenotypes on UK Biobank exome sequencing (or other) data. It uses [Snakemake](https://snakemake.github.io/) to parallelize large parts of the analysis across chromosomes and phenotypes. Conda on a Linux system is the main requirement for getting started.

This pipeline was originally run on different virtual machines in local mode. For full cluster support, users can easily modify `./run.sh` to use SLURM, SGE or other schedulers with Snakemake.

## Getting Started
### Installing basic dependencies
After cloning the repository, input files need to be configured and the basic installation script `./install.sh` needs to be run.

```
# download plink 1.9 and plink 2.0 , install snakemake in a separate local environment:
./install.sh
```
### Input file configuration
Paths to **phenotypes and covariates** are defined in `conf/config.yaml`.

Phenotypes and covariates are expected to be given in separate tab-delimeted files with single header lines that correspond to the names of variables. The first column in both files is expected to be called `iid` and contain the individual identifiers of participants. The files can be gzipped. These files are expected to be already pre-processed and filtered. Typically, this involves quantile transformation for the phenotypes, and the removal of related individuals. The specific entries that need to be modified are `covariates`, `phenotypes`, and `covariate_column_names`.

Paths to the **genotypes** in PLINK (.bim, .bam, .bed) format are defined in a separate file, defined with `samplesheet_plink` in config.yaml. We provide an example file at `conf/plinkfiles.tsv`.

> Note: The pipeline can't handle separate FIDs and IIDs, currently it is assumed that all FID == IID in your plink files (as is the case for UK Biobank data).

Some steps require a **reference genome and gene annotation**. We use [Ensembl](https://www.ensembl.org/index.html) gene annotations. Make sure all chromosome names are matched to your other data (e.g remove the "chr"-prefix where necessary). Adjust the following entries in the config file: `reference_fa` and `gene_annotation_gtf`. 

We used Ensembl 97 annotations available at ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/. If you want to use later versions, make sure to also adjust `ensembl_version` in config.yaml as well.

Another input to the pipeline which we can't provide ourselves are the [SpliceAI](https://github.com/Illumina/SpliceAI) variant effect predictions for SNVs. These can be downloaded from Illumina directly [here](https://basespace.illumina.com/s/otSPW8hnhaZR). Specifically, we used spliceai_scores.masked.snv.hg38.vcf.gz (v1.3). Make sure to adjust `splice_ai_scores_snv_vcf` in config.yaml.

### Running rules

Basic setup rules are defined in `./Snakemake`, whereas other rules are located in separate Snakemake modules within `./snake/`. It is recommended to execute rules using `./run.sh`. This script should be adjusted according to your cluster needs. *By default, everything is run locally*. To check if you have correctly configured your config file and initialize the work directory, run:

```
./run.sh setup_all
```

Users should adjust `./run.sh` if they wish for it to submit jobs to their cluster. Some of the rules in `./snake/deepripe.smk` require GPU support to perform predictions for RBP-binding.

The pipeline implements extensive variant annotation and association tests using functionally informed gene-based variant collapsing and kernel-based tests. As of March 2021, both the pipeline and the software it depends on to perform association tests ([seak](https://github.com/HealthML/seak)) are under active development. Forking is recommended.

### Variant effect prediction

The rules needed to run variant effect prediction are `evep_missense_proc_all` for missense variants, and `process_ensembl_vep_output_highimpact_all` for protein loss of function variants (pLOF). The rule `splice_ai_filter_and_overlap_with_genotypes_all` overlaps all SNVs to published SpliceAI predictions. Finally, `run_deepripe_vep_all` runs DeepRiPe variant effect predictions for the binding of RBPs.

### Association tests

To run a basic analysis using gene-based variant collapsing score tests for missense and pLOF variants, use the rule `assoc_baseline_scoretest_all`. This will also export (tiny!) HDF5-files containing the pre-processed burdens for all individuals which have complete covariates to `work/gene_scores/indicator01/missense/all/` and `work/gene_scores/indicator01/pLOF/all/`.

To perform kernel-based tets with the sLRT, the follwing rules can be used: To run gene-based variant collapsing and kernel-based tests for for missense variants use rule `assoc_missense_localcollapsing_all`. To run gene-based variant collapsing and kernel-based tests for for splice-variants use rule `assoc_spliceai_linw_all`. To run kernel-based tests incorporating variant effect predictions for RBP-binding use rule `assoc_deepripe_multiple_cholesky_all`.
