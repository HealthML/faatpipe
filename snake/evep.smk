
rule install_ensembl_cache:
    # use this rule to install a specific version of the ensembl VEP cache. Ideally this version should match that of the gene annotation GTF.
    output:
        expand('data/evep_cache/homo_sapiens/{ensembl_version}_GRCh38/info.txt', ensembl_version=config['ensembl_version'])
    params:
        ensembl_version = config['ensembl_version']
    conda:
        '../env/ensembl_vep.yml'
    log:
        expand("logs/evep/install_cache_{ensembl_version}.log", ensembl_version=config['ensembl_version'])
    shell:
        "("
        "rm -r $(dirname {output}); "
        "vep_install -a cf -s homo_sapiens -y GRCh38 -c ./data/evep_cache/ --CONVERT --CACHE_VERSION {params.ensembl_version}"
        ") &> {log} "
        

rule run_ensembl_vep:
    input:
        cache_info=expand('data/evep_cache/homo_sapiens/{ensembl_version}_GRCh38/info.txt', ensembl_version=config['ensembl_version']),
        vcf = rules.bim_to_vcf.output.vcf
    output:
        vep_tsv = 'work/variant_effect_prediction/ensembl_vep/{id}.tsv.gz',
        vep_html = 'work/variant_effect_prediction/ensembl_vep/{id}.html'
    params:
        ensembl_version=config['ensembl_version']
    conda:
        '../env/ensembl_vep.yml'
    log:
        'logs/evep/run_ensembl_vep_{id}.log'
    shell:
        "("
        "vep --format vcf "
        "--offline "
        "--cache "
        "--dir_cache data/evep_cache/ "
        "--assembly GRCh38 "
        "--species homo_sapiens "
        "--cache_version {params.ensembl_version} "
        "--input_file {input.vcf} "
        "--sift b "
        "--polyphen b "
        "--output_file {output.vep_tsv} "
        "--stats_file {output.vep_html} "
        "--compress_output gzip "
        ") &> {log} "
        # --fork 4
        # we don't use forking because we are already parallelizing over chromosomes
        
        
rule run_ensembl_vep_all:
    # rule that triggers running the ensembl variant effect predictor for all chromosomes 
    input:
        expand(rules.run_ensembl_vep.output.vep_tsv, id=plinkfiles.getIds())
    output:
        touch('work/variant_effect_prediction/ensembl_vep/all.ok')
        
        
rule filter_ensembl_vep_output_highimpact:
    # get's all the variants marked with IMPACT=HIGH
    input:
        tsv = rules.run_ensembl_vep.output.vep_tsv
    output:
        tsv = 'work/variant_effect_prediction/ensembl_vep/filtered/high_impact/{id}.tsv.gz'
    log:
        "logs/evep/filter_ensembl_vep_output_highimpact_{id}.log"
    shell:
        "("
        "bash script/bash/filter_ensembl_vep_output_highimpact.sh {input.tsv}  | gzip > {output.tsv} "
        ") &> {log}"


rule process_ensembl_vep_output_highimpact:
    # filters the output of filter_ensembl_vep_output_highimpact so there's only one entry for every variant
    input:
        tsv = rules.filter_ensembl_vep_output_highimpact.output.tsv
    output:
        tsv = 'work/variant_effect_prediction/ensembl_vep/processed/high_impact/{id}.tsv.gz'
    log: 
        'logs/evep/process_ensembl_vep_output_highimpact_{id}.log'
    script:
        '../script/python/evep_process_high_impact_variants.py'
        

rule process_ensembl_vep_output_highimpact_all:
    # runs high-impact (LOF) variant effect prediction for all chromosomes
    input:
        expand(rules.process_ensembl_vep_output_highimpact.output, id=plinkfiles.getIds())


rule filter_ensembl_vep_output_missense:
    # get's all the missense variants
    # this is a dead end, not actually needed!
    input:
        tsv = rules.run_ensembl_vep.output.vep_tsv
    output:
        tsv = 'work/variant_effect_prediction/ensembl_vep/filtered/missense/raw/{id}.tsv.gz'
    log:
        "logs/evep/filter_ensembl_vep_output_missense_{id}.log"
    shell:
        "("
        "bash script/bash/filter_ensembl_vep_output_missense.sh {input.tsv}  | gzip > {output.tsv} "
        ") &> {log}"
        
        
rule download_ensembl_cds_fasta:
    # downloads the CDS reference fasta so we can get the AA-trimers of missense mutations
    # the awk one-liner clips away unnecessary suffixes of transcript names
    params:
        download_link = 'ftp://ftp.ensembl.org/pub/release-{}/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz'.format(config['ensembl_version'])
    output:
        fasta = 'data/reference/cds/Homo_sapiens.GRCh38.cds.all.rename.fa',
        fai = 'data/reference/cds/Homo_sapiens.GRCh38.cds.all.rename.fa.fai'
    log:
        'logs/evep/download_ensembl_cds_fasta.log'
    conda:
        '../env/genomics.yml'
    shell:
        "("
        "wget --directory-prefix=data/reference/cds/ {params.download_link} ; "
        "zcat data/reference/cds/Homo_sapiens.GRCh38.cds.all.fa.gz | awk '{{if($1 ~ />/){{print substr($1,1,16)}}else{{print $0}}}}' > {output.fasta}; "
        "rm data/reference/cds/Homo_sapiens.GRCh38.cds.all.fa.gz ; "
        "samtools faidx {output.fasta} "
        ") &> {log}"
        
    
rule evep_missense_proc:
    # rule to run missense variant annotation for a single chromosome
    input:
        tsv = rules.run_ensembl_vep.output.vep_tsv,
        cds_fasta = rules.download_ensembl_cds_fasta.output.fasta
    output:
        tsv_all = 'work/variant_effect_prediction/ensembl_vep/processed/missense/{id}.tsv.gz',
        tsv_filtered = 'work/variant_effect_prediction/ensembl_vep/processed/missense/{id}_filtered.tsv'
    log:
        'logs/evep/evep_missense_proc_{id}.log'
    conda:
        '../env/genomics.yml'
    shell:
        "("
        "python script/python/evep_missense_get_3grams.py "
        "--infile {input.tsv} "
        "--out_prefix work/variant_effect_prediction/ensembl_vep/processed/missense/{wildcards.id} "
        "--cds_fasta {input.cds_fasta} && gzip work/variant_effect_prediction/ensembl_vep/processed/missense/{wildcards.id}.tsv ;"
        ") &> {log}"
    
        
rule evep_missense_proc_all:
    # rule to run missense variant annotation and filtering for all chromosomes
    input:
        expand(rules.evep_missense_proc.output, id=plinkfiles.getIds())



rule export_plof_burden:
    # export protein LOF indicator variables to (tiny!) files
    input:
        mac_report = rules.mac_report.output.tsv,
        anno_tsv = rules.basic_annotation.output.tsv,
        complete_cases = 'data/covariates/complete_cases/covariates.txt',
        vep_tsv = rules.process_ensembl_vep_output_highimpact.output.tsv,
        genotypes_bed = rules.link_genotypes.output.bed,
        regions_bed = rules.get_protein_coding_genes.output.pc_genes_bed,
        install_seak = rules.install_seak.output
    params:
        effect = 'LOF',
        max_missing = config['max_missingness'],
        max_maf = config['maf_cutoff'],
        filter_highconfidence = lambda wc: {'all':False, 'highconf_only':True}[wc.filter_highconfidence]
    output:
        h5 = 'work/gene_scores/indicator01/pLOF/{filter_highconfidence}/{id}.h5',
        gene_txt = 'work/gene_scores/indicator01/pLOF/{filter_highconfidence}/{id}_gene.txt',
        iid_txt = 'work/gene_scores/indicator01/pLOF/{filter_highconfidence}/{id}_iid.txt'
    log:
        'logs/evep/export_plof_burden_{filter_highconfidence}_{id}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/export_burden.py'
        
        
rule export_plof_burden_all:
    # run the rule above for all chromosomes
    input:
        expand(rules.export_plof_burden.output, id = plinkfiles.getIds(), filter_highconfidence=['all']) # filter_highconfidence=['all','highconf_only'])


rule export_missense_burden:
    # export missense indicator variables to (tiny!) files
    # Note: these are *not* the ones used in association tests in the publication!
    input:
        mac_report = rules.mac_report.output.tsv,
        anno_tsv = rules.basic_annotation.output.tsv,
        complete_cases = 'data/covariates/complete_cases/covariates.txt',
        vep_tsv = rules.evep_missense_proc.output.tsv_filtered,
        genotypes_bed = rules.link_genotypes.output.bed,
        regions_bed = rules.get_protein_coding_genes.output.pc_genes_bed,
        install_seak = rules.install_seak.output
    params:
        effect = 'missense',
        max_missing = config['max_missingness'],
        max_maf = config['maf_cutoff'],
        min_impact = config['min_impact'],
        filter_highconfidence = lambda wc: {'all':False, 'highconf_only':True}[wc.filter_highconfidence]
    output:
        h5 = 'work/gene_scores/indicator01/missense/{filter_highconfidence}/{id}.h5',
        gene_txt = 'work/gene_scores/indicator01/missense/{filter_highconfidence}/{id}_gene.txt',
        iid_txt = 'work/gene_scores/indicator01/missense/{filter_highconfidence}/{id}_iid.txt'
    log:
        'logs/evep/export_missense_burden_{filter_highconfidence}_{id}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/export_burden.py'

rule export_missense_burden_all:
    # run the rule above for all chromosomes
    input:
        expand(rules.export_missense_burden.output, id = plinkfiles.getIds(), filter_highconfidence=['all'])
        
