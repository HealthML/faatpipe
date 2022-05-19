
from snake.utils import SampleSheet
from script.python.util.snake import clean_str

import os
import gzip

configfile: "conf/config.yaml"

###################################
# samplesheet for the plink files #
###################################

# load plink files
plinkfiles = SampleSheet(config['samplesheet_plink'])

for col in ['bim','fam','bed']:
    assert col in plinkfiles.sheet.columns, 'Error: column {} is missing from samplesheet!'.format(col)
    
    
##############
# phenotypes #
##############

if 'phenotypes' in config:
    if config['phenotypes'].endswith('.gz'):
        with gzip.open(config['phenotypes'],'rt') as f:
            phenotypes = f.readline().rstrip().split('\t')
    else:
        with open(config['phenotypes'], 'r') as f:
            phenotypes = f.readline.rstrip().split('\t')
    # "phenotypes" provides a mapping from file-prefixes to columns in the phenotype.tsv file
    phenotypes = { clean_str(p): p  for p in phenotypes[1:] }


########################
# wildcard constraints #
########################

wildcard_constraints:
    id="[a-z]*[0-9]+",
    strand="plus|minus"


###############
# local rules #
###############

localrules: link_all, install_sh, install_seak, setup_all, run_deepripe_vep_all, run_ensembl_vep_all, evep_missense_proc_all, splice_ai_filter_and_overlap_with_genotypes_all, splice_ai_vcf_to_tsv_all, filter_variants_all, export_plof_burden_all, export_missense_burden_all, assoc_baseline_scoretest_all, mac_report_all, assoc_spliceai_linw_all, assoc_deepripe_single_localcollapsing_all, pLOF_nsnp_cummac_all


###############
# dummy rule  #
###############

rule dummy:
    # dummy rule to test things
    output:
        temp(touch("dummy.ok"))
    conda:
        'env/ensembl_vep.yml'
        #'env/seak.yml'
    shell:
        "echo success."


################
#  setup rules #
################

rule link_genotypes:
    input:
        bim = lambda wc: plinkfiles.get('bim', wc.id),
        fam = lambda wc: plinkfiles.get('fam', wc.id),
        bed = lambda wc: plinkfiles.get('bed', wc.id)
    output:
        bim = "data/genotypes/{id}.bim",
        fam = "data/genotypes/{id}.fam",
        bed = "data/genotypes/{id}.bed"
    run:
        if not os.path.isdir('data/genotypes'):
            os.makedirs('data/genotypes')
        shell('ln -s -r "{input.bim}" {output.bim} && ln -s -r "{input.fam}" {output.fam} &&  ln -s -r "{input.bed}" {output.bed}')

rule link_reference:
    input:
        config['reference_fa']
    output:
        'data/reference/genome.fa'
    run:
        if not os.path.isdir('data/reference'):
            os.makedirs('data/reference')
        shell('ln -s -r "{input}" {output}')

rule link_gene_annotation:
    input:
        config['gene_annotation_gtf']
    output:
        'data/reference/gene_annotation.gtf' + '.gz' if config['gene_annotation_gtf'].endswith('.gz') else ''
    run:
        if not os.path.isdir('data/reference'):
            os.makedirs('data/reference')
        shell('ln -s -r "{input}" {output}')

rule link_all:
    # checks if all the input files are present and links some of them to the working directory
    input:
        rules.link_gene_annotation.output,
        rules.link_reference.output,
        expand('data/genotypes/{id}.{ext}', id=plinkfiles.getIds(), ext=['bim','fam','bed']),
        config['high_confidence_regions'],
        config['exome_sequencing_target_regions'],
        #config['iid_withdraw'],
        #config['iid_exclude']
    output:
        touch('data/workdir_init.ok')
        
        
rule install_sh:
    # runs install.sh
    # this could do a lot of stuff, but right now it just downloads plink1.9, plink2 and installs snakemake in a separate environment used by run.sh
    output:
        plink2 = 'bin/plink2',
        plink1 = 'bin/plink'
    shell:
        './install.sh'

rule install_seak:
    # pulls the specific commit that was used at time of publication.
    output:
        touch('install/install_seak.ok')
    conda:
        'env/seak.yml'
    shell:
        'if [ ! -d seak ]; then '
        'git clone https://github.com/HealthML/seak.git && cd seak && git checkout tags/v0.4.2 ; '
        'else '
        'cd seak && git checkout tags/v0.4.2 ; '
        'fi ; '
        'pip install -e . '
        

rule setup_all:
    input:
        rules.link_all.output,
        rules.install_sh.output
    output:
        touch('data/install.ok')
    
#####################
# snakemake modules #
#####################

include: "snake/setup.smk"
include: "snake/deepripe.smk"
include: "snake/evep.smk"
include: "snake/splice_ai.smk"
include: "snake/association.smk"
include: "snake/results_processing.smk"
