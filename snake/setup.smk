

# These rules perform a couple of format conversions, and generate basic statistics on the variants analysed
# the {id} wildcard refers to the different chromosomes
# the {pheno} wildcard refers to the different phenotypes

rule get_protein_coding_genes:
    # get's the protein coding genes from the annotation GTF file 
    input:
        gtf = config['gene_annotation_gtf'],
        target_regions = config['exome_sequencing_target_regions']
    output:
        pc_genes_bed = 'data/reference/gene_annotation/protein_coding_genes.bed',
        pc_genes_bed_plus = 'data/reference/gene_annotation/protein_coding_genes_plus.bed',
        pc_genes_bed_minus = 'data/reference/gene_annotation/protein_coding_genes_minus.bed'
    log:
        'logs/setup/get_protein_coding_genes.log'
    conda:
        '../env/genomics.yml'
    script:
        '../script/python/get_protein_coding_genes_from_gtf.py'


rule bim_to_vcf:
     # converts the bim file(s) to VCF format, which is needed for most variant effect prediction steps
     # there should be a better way to do this but the plink documentation is so bad I can't figure it out.
     input:
         bim = 'data/genotypes/{id}.bim',
         fam = 'data/genotypes/{id}.fam',
         plink2 = 'bin/plink2'
     output:
         vcf = 'work/variant_effect_prediction/vcf/{id}.vcf.gz'
     params:
         in_prefix = lambda wc, input: input.bim.replace('.bim',''),
         out_prefix = lambda wc, output: output.vcf.replace('.vcf.gz','')
     log:
         'logs/setup/bim_to_vcf_{id}.log'
     shell:
         "("
         "head -n 1 {input.fam} | cut -f1 > {params.out_prefix}_tmp_fam.txt; "
         "{input.plink2} --memory 1000 --bfile {params.in_prefix} --export vcf-4.2 id-paste=iid bgz --keep-fam {params.out_prefix}_tmp_fam.txt "
         "--out {params.out_prefix} && rm {params.out_prefix}_tmp_fam.txt "
         ") &> {log}"
         

rule bim_to_vcf_all:
    # rule to run bim -> vcf for all samples
    input:
        expand(rules.bim_to_vcf.output.vcf, id=plinkfiles.getIds())
        
                
rule basic_annotation:
    # overlap variants with gene annotation and
    # overlap variants with high-confidence regions and
    # overlap variants with exome-sequencing target regions and
    # export overlaps to TSV
    input:
        vcf = rules.bim_to_vcf.output.vcf,
        hc_bed = config['high_confidence_regions'],
        es_bed = config['exome_sequencing_target_regions'],
        anno_gtf = rules.link_gene_annotation.output[0]
    output:
        tsv = 'work/basic_annotation/{id}.tsv.gz'
    log:
        'logs/setup/basic_annotation_{id}.log'
    conda:
        '../env/genomics.yml'
    script:
        '../script/python/basic_annotation.py'
        
    
rule gather_complete_cases:
    # creates txt files containing the complete cases for the phenotypes / covariates.
    input:
        covar_tsv = config['covariates'],
        pheno_tsv = config['phenotypes']
    output:
        expand('data/covariates/complete_cases/{pheno}.txt', pheno=list(phenotypes.keys()) + ['covariates'])
    log:
        'logs/setup/gather_complete_cases.log'
    script:
        '../script/python/gather_complete_cases.py'
        


rule complete_cases_ancestry:
    # create ancestry keep-files for plink
    input:
        iid_txt = 'data/covariates/complete_cases/covariates.txt',
        ancestry_tsv = config['ancestry_scoring_file']
    output:
        ancestry_keep = expand('data/ancestry_keep_files/{ancestries}.keep', ancestries=['AFR','AMR','EAS','EUR','SAS'])
    run:
        import pandas as pd
        
        with open(input['iid_txt'], 'r') as infile:
            iids = [ l.rstrip() for l in infile ]
            
        ancestry = pd.read_csv(input['ancestry_tsv'], sep='\t', dtype={0:str})
        if 'IID' in ancestry.columns:
            ancestry.rename(columns={'IID':'iid'}, inplace=True)
        ancestry.set_index('iid', inplace=True)    
        
        # print(ancestry.head())
        
        print('{} individuals in ancestry file'.format(len(ancestry)))
        print('{} individuals in iid-file. Subsetting to those individuals.'.format(len(iids)))
        
        ancestry = ancestry.loc[iids]

        for a in ancestry.columns:
            outfile = 'data/ancestry_keep_files/{}.keep'.format(a)
            
            with open(outfile, 'w') as out:
                for idx in ancestry[ancestry[a] > 0.5].index.values:
                    out.write('{}\t{}\n'.format(idx, idx))
            
            print('written {} iids to {}'.format((ancestry[a] > 0.5).sum(), outfile))
            
    
rule mac_report:
    # filter by HWE and
    # count minor allele frequencies for covariate-complete cases
    input:
        # anno_tsv = rules.basic_annotation.output.tsv,
        iid_txt = 'data/covariates/complete_cases/covariates.txt',
        plink = 'bin/plink',
        bed = rules.link_genotypes.output.bed
    output:
        frq_tsv = temp('work/mac_report/all/{id}.frq.counts.gz'),
        tsv = 'work/mac_report/all/{id}.tsv.gz',
        plink_log = 'work/mac_report/all/{id}.log'
    params:
        in_prefix = lambda wc, input: input.bed[:-4],
        out_prefix = 'work/mac_report/all/{id}'
    log:
        'logs/mac_report/all/{id}.log'
    conda:
        '../env/genomics.yml'
    threads:
        1
    shell:
        'python script/python/mac_report.py '
        '--iid {input.iid_txt} '
        '--bed {params.in_prefix} '
        '-o {params.out_prefix} '
        '--plink_path {input.plink} '
        '--log {log} '
        '--threads {threads} '
        '&> {log} '

rule mac_report_all:
    input:
        expand(rules.mac_report.output, id = plinkfiles.getIds())
        
        
rule mac_report_ancestry:
    # filter by HWE and
    # create minor allele count report for specific ancestry and chromosome
    input:
        iid_txt = 'data/ancestry_keep_files/{ancestry}.keep',
        plink = 'bin/plink',
        bed = rules.link_genotypes.output.bed
    output:
        frq_tsv = temp('work/mac_report/all/{ancestry}/{id}.frq.counts.gz'),
        tsv = 'work/mac_report/all/{ancestry}/{id}.tsv.gz',
        plink_log = 'work/mac_report/all/{ancestry}/{id}.log'
    params:
        in_prefix = lambda wc, input: input.bed[:-4],
        out_prefix = 'work/mac_report/all/{ancestry}/{id}'
    log:
        'logs/mac_report_ancestry/{ancestry}/{id}.log'
    conda:
        '../env/genomics.yml'
    threads:
        1
    shell:
        'python script/python/mac_report.py '
        '--iid {input.iid_txt} '
        '--bed {params.in_prefix} '
        '-o {params.out_prefix} '
        '--plink_path {input.plink} '
        '--log {log} '
        '--threads {threads} '
        '&> {log} '
        

rule mac_report_ancestry_all:
    # create minor allele count reports for all ancestries and all chromosomes
    input:
        expand(rules.mac_report_ancestry.output, id = plinkfiles.getIds(), ancestry=['AFR','AMR','EAS','EUR','SAS'])

        
rule filter_variants:
    # apply missingness filters and
    # high-confidence region filter for phenotype-complete cases
    # returns the filtered variant ids, their MAF (calculated on the covariate-complete cases) and their MAC (calculated on the *phenotype*-complete cases)
    # the MAF filter is applied while performing the association tests, NOT here.
    input:
        mac_tsv = rules.mac_report.output.tsv,
        anno_tsv = rules.basic_annotation.output.tsv,
        iid_txt = 'data/covariates/complete_cases/{pheno}.txt',
        bed = rules.link_genotypes.output.bed,
        plink = 'bin/plink',
    output:
        vid_tsv = 'work/mac_report/{pheno}/{id}.tsv.gz'
    params:
        in_prefix = lambda wc, input: input.bed[:-4],
        out_prefix = 'work/mac_report/{pheno}/{id}',
        max_missing = config['max_missingness']
    log:
        'logs/mac_report/{pheno}/{id}.log'
    conda:
        '../env/genomics.yml'
    threads:
        1
    shell:
        '('
        'python script/python/filter_variants.py '
        '--bed {params.in_prefix} '
        '--iid {input.iid_txt} '
        '--mac_report {input.mac_tsv} '
        '--anno {input.anno_tsv} '
        '--max_missing {params.max_missing} '
        '--plink_path {input.plink} '
        '--out_prefix {params.out_prefix} '
        '--log {log} '
        '--threads {threads} ' 
        ') &> {log} '
        
        
rule filter_variants_all:
    # runs the rule above for all chromosomes and phenotypes
    input:
        expand(rules.filter_variants.output.vid_tsv, id = plinkfiles.getIds(), pheno=phenotypes.keys())


