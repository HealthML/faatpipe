
# TODO: add conditional_list to config.yaml


rule assoc_baseline_scoretest_conditional_analysis:
    # performs conditional analysis for the baseline 
    input:
        h5_missense = rules.export_missense_burden.output.h5,
        iid_missense = rules.export_missense_burden.output.iid_txt,
        gid_missense = rules.export_missense_burden.output.gene_txt,
        h5_lof = rules.export_plof_burden.output.h5,
        iid_lof = rules.export_plof_burden.output.iid_txt,
        gid_lof = rules.export_plof_burden.output.gene_txt,
        covariates_tsv = config['covariates'],
        phenotypes_tsv = config['phenotypes'],
        results_tsv = rules.assoc_baseline_scoretest.output['results_tsv'],
        conditional_list = config['conditional_list'],
        conditional_geno = config['conditional_geno'] + '.bed'
    output:
        results_tsv = 'work/association/baseline_scoretest/{filter_highconfidence}/{pheno}/conditional_analysis/conditionalres_{id}.tsv.gz'
    params:
        # phenotypes is a dictionary that maps file-prefixes to column names, see ../Snakefile
        phenotype = lambda wc: phenotypes[ wc.pheno ],
        covariate_column_names = config['covariate_column_names'],
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        columns = ['pv_pLOF','pv_missense','pv_mrg'],
        significance_cutoff=1e-7
    threads:
        1
    resources:
        mem_mb=4000,
        time='0:10:00'
    log:
        'logs/association/baseline_scoretest_conditional_analysis/{filter_highconfidence}/{pheno}/{id}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_baseline_scoretest_conditional_analysis.py'
        
        
rule assoc_baseline_scoretest_conditional_analysis_all:
    input:
        expand(rules.assoc_baseline_scoretest_conditional_analysis.output, filter_highconfidence=['all'], pheno=phenotypes.keys(), id=plinkfiles.getIds())
        

rule assoc_missense_localcollapsing_conditional_analysis:
    # calculates gene-specific LRT p-values for the most significant genes for missense variants
    # conditions on a list of common variants
    input:
        h5_lof = expand(rules.export_plof_burden.output.h5, id = plinkfiles.getIds(), allow_missing=True),
        iid_lof = expand(rules.export_plof_burden.output.iid_txt,  id = plinkfiles.getIds(), allow_missing=True),
        gid_lof = expand(rules.export_plof_burden.output.gene_txt, id = plinkfiles.getIds(), allow_missing=True),
        covariates_tsv = config['covariates'],
        phenotypes_tsv = config['phenotypes'],
        bed = expand(rules.link_genotypes.output.bed, id = plinkfiles.getIds()),
        vep_tsv = expand(rules.evep_missense_proc.output.tsv_filtered, id = plinkfiles.getIds()),
        mac_report = expand(rules.filter_variants.output.vid_tsv, id = plinkfiles.getIds(), allow_missing=True),
        regions_bed = rules.get_protein_coding_genes.output.pc_genes_bed,
        results_tsv = rules.assoc_missense_localcollapsing.output.results_tsv,
        seak_install = rules.install_seak.output,
        conditional_list = config['conditional_list'],
        conditional_geno = config['conditional_geno'] + '.bed'
    output:
        results_tsv = 'work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/conditional_analysis/lrt_conditional.tsv.gz'
    params:
        phenotype = lambda wc: phenotypes[ wc.pheno ],
        covariate_column_names = config['covariate_column_names'],
        kernels = ['linwcollapsed','linwcollapsed_cLOF','linwb','linwb_mrgLOF'], # kernels to consider
        max_maf = config['maf_cutoff'],
        min_impact = config['min_impact'],
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        out_dir_stats=lambda wc, output: '/'.join(output['results_tsv'].split('/')[:-1]),
        significance_cutoff=1e-7,
        debug=False
    resources:
        mem_mb=12000,
        time="02:30:00"
    log:
        'logs/association/sclrt_kernels_missense_conditional_analysis/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_missense_conditional_analysis.py'


rule all_assoc_missense_localcollapsing_conditional_analysis:
    # runs rule above for all phenotypes
    input:
        expand(rules.assoc_missense_localcollapsing_conditional_analysis.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])
        
        
rule assoc_spliceai_linw_conditional_analysis:
    # calculates gene-specific LRT p-values for the most significant genes for missense variants
    # conditions on a list of common variants
    input:
        h5_lof = expand(rules.export_plof_burden.output.h5, id = plinkfiles.getIds(), allow_missing=True),
        iid_lof = expand(rules.export_plof_burden.output.iid_txt,  id = plinkfiles.getIds(), allow_missing=True),
        gid_lof = expand(rules.export_plof_burden.output.gene_txt, id = plinkfiles.getIds(), allow_missing=True),
        ensembl_vep_tsv = expand(rules.process_ensembl_vep_output_highimpact.output.tsv, id = plinkfiles.getIds(), allow_missing=True),
        covariates_tsv = config['covariates'],
        phenotypes_tsv = config['phenotypes'],
        bed = expand(rules.link_genotypes.output.bed, id = plinkfiles.getIds()),
        vep_tsv = expand(rules.splice_ai_filter_and_overlap_with_genotypes.output.tsv, id = plinkfiles.getIds()),
        mac_report = expand(rules.filter_variants.output.vid_tsv, id = plinkfiles.getIds(), allow_missing=True),
        regions_bed = rules.get_protein_coding_genes.output.pc_genes_bed,
        results_tsv = rules.assoc_spliceai_linw.output.results_tsv,
        seak_install = rules.install_seak.output,
        conditional_list = config['conditional_list'],
        conditional_geno = config['conditional_geno'] + '.bed'
    output:
        results_tsv = 'work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/conditional_analysis/lrt_conditional.tsv.gz'
    params:
        kernels = ['linwb','linw','linwb_mrgLOF','linw_cLOF'],
        phenotype = lambda wc: phenotypes[wc.pheno],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['splice_ai_min_impact'],
        out_dir_stats=lambda wc, output: '/'.join(output['results_tsv'].split('/')[:-1]),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        debug = False,
        significance_cutoff = 1e-7
    resources:
        mem_mb=10000,
        time="02:30:00"
    log:
        'logs/association/sclrt_kernels_spliceai_conditional_analysis/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_spliceai_conditional_analysis.py'



rule all_assoc_spliceai_linw_conditional_analysis:
    # runs rule above for all phenotypes
    input:
        expand(rules.assoc_spliceai_linw_conditional_analysis.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])



rule assoc_deepripe_multiple_cholesky_conditional_analysis:
    # re-tests most significant genes for DeepRiPe variants using gene-specific null-distributions
    input:
        bed = expand(rules.link_genotypes.output.bed, id = plinkfiles.getIds()),
        h5_rbp_plus = expand(rules.run_deepripe_vep.output.h5, id = plinkfiles.getIds(), strand=['plus']),
        h5_rbp_minus = expand(rules.run_deepripe_vep.output.h5, id = plinkfiles.getIds(), strand=['minus']),
        bed_rbp_plus = expand(rules.run_deepripe_vep.output.bed, id = plinkfiles.getIds(), strand=['plus']),
        bed_rbp_minus = expand(rules.run_deepripe_vep.output.bed, id = plinkfiles.getIds(), strand=['minus']),
        ensembl_vep_tsv = expand(rules.process_ensembl_vep_output_highimpact.output.tsv, id = plinkfiles.getIds(), allow_missing=True),
        mac_report = expand(rules.filter_variants.output.vid_tsv, id = plinkfiles.getIds(), allow_missing=True),
        regions_bed = rules.get_protein_coding_genes.output.pc_genes_bed,
        seak_install = rules.install_seak.output,
        covariates_tsv = config['covariates'],
        phenotypes_tsv = config['phenotypes'],
        results_tsv = rules.assoc_deepripe_multiple_cholesky.output.results_tsv,
        conditional_list = config['conditional_list'],
        conditional_geno = config['conditional_geno'] + '.bed'
    output:
        results_tsv = 'work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/conditional_analysis/lrt_conditional.tsv.gz'
    params:
        kernels = ['linwcholesky', 'linwcholesky_notLOF'],
        phenotype = lambda wc: phenotypes[wc.pheno],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['deepripe_min_impact'],
        out_dir_stats=lambda wc, output: '/'.join(output['results_tsv'].split('/')[:-1]),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        rbp_of_interest = rules.assoc_deepripe_multiple_cholesky.params.rbp_of_interest,
        debug = False,
        random = False,
        significance_cutoff = 1e-7
    resources:
        mem_mb=10000,
        time="02:30:00"
    log:
        'logs/association/sclrt_kernels_deepripe_multiple_conditional_analysis/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_deepripe_multiple_conditional_analysis.py'


rule assoc_deepripe_multiple_cholesky_conditional_analysis_all:
    # runs rule above for all phenotypes - this will effectively run the entire pipeline for DeepRiPe variants
    input:
        expand(rules.assoc_deepripe_multiple_cholesky_conditional_analysis.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])
        
        
rule conditional_analysis_all:
    input:
        rules.assoc_deepripe_multiple_cholesky_conditional_analysis_all.input,
        rules.all_assoc_spliceai_linw_conditional_analysis.input,
        rules.all_assoc_missense_localcollapsing_conditional_analysis.input,
        rules.assoc_baseline_scoretest_conditional_analysis_all.input