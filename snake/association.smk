
# These rules implement association tests 
# the {id} wildcard refers to the different chromosomes
# the {pheno} wildcard refers to the different phenotypes


rule assoc_baseline_scoretest:
    # runs "baseline" association tests:
    # score tests for indicator variables
    # 1 if at least one loss of function variant present, 0 otherwise
    # 1 if at least one high-impact missense variant present, 0 otherise 
    # "high-impact" for missense variants is defined buy the "min_impact" parameter in the config.
    input:
        h5_missense = rules.export_missense_burden.output.h5,
        iid_missense = rules.export_missense_burden.output.iid_txt,
        gid_missense = rules.export_missense_burden.output.gene_txt,
        h5_lof = rules.export_plof_burden.output.h5,
        iid_lof = rules.export_plof_burden.output.iid_txt,
        gid_lof = rules.export_plof_burden.output.gene_txt,
        covariates_tsv = config['covariates'],
        phenotypes_tsv = config['phenotypes']
    params:
        # phenotypes is a dictionary that maps file-prefixes to column names, see ../Snakefile
        phenotype = lambda wc: phenotypes[ wc.pheno ],
        covariate_column_names = config['covariate_column_names'],
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence]
    output:
        results_tsv = 'work/association/baseline_scoretest/{filter_highconfidence}/{pheno}/results_{id}.tsv.gz'
    threads:
        1
    resources:
        mem_mb=4000,
        time='2:00:00'
    log:
        'logs/association/baseline_scoretest/{filter_highconfidence}/{pheno}/{id}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_baseline_scoretest.py'
                
        
rule assoc_baseline_scoretest_all:
    input:
        expand(rules.assoc_baseline_scoretest.output, pheno=phenotypes.keys(), id=plinkfiles.getIds(), filter_highconfidence=['all'])
      
      
rule pLOF_nsnp_cummac:
    input:
        basic_anno_tsv=expand(rules.basic_annotation.output.tsv, id=plinkfiles.getIds()),
        mac_report=expand(rules.filter_variants.output.vid_tsv, id=plinkfiles.getIds(), allow_missing=True),
        evep_tsv = expand(rules.process_ensembl_vep_output_highimpact.output.tsv, id=plinkfiles.getIds()),
        pc_genes = rules.get_protein_coding_genes.output.pc_genes_bed
    output:
        tsv='work/association/baseline_scoretest/{filter_highconfidence}/{pheno}/pLOF_nSNP_cumMAC.tsv.gz'
    conda:
        '../env/seak.yml'
    resources:
        mem_mb=4000,
        time='0:30:00'
    script:
        '../script/python/pLOF_nSNP_cumMAC.py'
      
rule pLOF_nsnp_cummac_all:
    input:
        expand(rules.pLOF_nsnp_cummac.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])
        

rule assoc_missense_localcollapsing:
    # runs association tests for missense variants using local collapsing
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
        seak_install = rules.install_seak.output
    output:
        results_tsv = 'work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/results.tsv.gz',
        chi2param_tsv = 'work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/chi2mixparam.tsv.gz'
    params:
        phenotype = lambda wc: phenotypes[ wc.pheno ],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/lrt_stats/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        sclrt_nominal_significance_cutoff = 0.1,
        debug=False
    log:
        'logs/association/sclrt_kernels_missense/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_missense.py'

rule assoc_missense_localcollapsing_all:
    # runs rule above for all phenotypes
    input:
        expand(rules.assoc_missense_localcollapsing.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])

rule assoc_missense_localcollapsing_eval_top_hits:
    # calculates single-variant p-values and regression coefficients and other statistics for the most significant genes for missense variants
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
        seak_install = rules.install_seak.output
    output:
        out_ok = touch('work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/top_hits/all.ok')
    params:
        kernels = ['lincollapsed','linwcollapsed','lincollapsed_cLOF','linwcollapsed_cLOF','linb','linwb','linb_mrgLOF','linwb_mrgLOF'], # genes with p < 1e-7 in one of these kernels will be analysed in detail
        phenotype = lambda wc: phenotypes[ wc.pheno ],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/top_hits/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence]
    log:
        'logs/association/sclrt_kernels_missense_eval_top_hits/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_missense_eval_top_hits.py'
        
        
rule assoc_missense_localcollapsing_eval_top_hits_all:
    # run above rule for all phenotypes - this rule will effectively run the entire pipeline for missense variants
    input:
        expand(rules.assoc_missense_localcollapsing_eval_top_hits.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])
        

rule assoc_missense_localcollapsing_retest_top_hits:
    # calculates gene-specific LRT p-values for the most significant genes for missense variants
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
        seak_install = rules.install_seak.output
    output:
        out_ok = touch('work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/lrtsim/all.ok'),
        results_tsv = 'work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/lrtsim/lrt_retest.tsv.gz'
    params:
        significance_cutoff = 5e-6, # genes with any p-value below this threshold will be analysed in detail...
        kernels = ['linwcollapsed','linwcollapsed_cLOF','linwb','linwb_mrgLOF'], # kernels to consider 
        phenotype = lambda wc: phenotypes[ wc.pheno ],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/lrtsim/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        debug=False,
        random=False
    log:
        'logs/association/sclrt_kernels_missense_retest_top_hits/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_missense_retest_top_hits.py'
        
        
rule assoc_missense_localcollapsing_retest_top_hits_all:
    # run above rule for all phenotypes
    input:
        expand(rules.assoc_missense_localcollapsing_retest_top_hits.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])
        

rule assoc_missense_localcollapsing_retest_random:
    # calculates gene-specific LRT p-values for random genes for missense variants
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
        seak_install = rules.install_seak.output
    output:
        out_ok = touch('work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/lrtsim_random/all.ok'),
        results_tsv = 'work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/lrtsim_random/lrt_retest.tsv.gz'
    params:
        significance_cutoff = 5e-6, # genes with any p-value below this threshold will be analysed in detail...
        kernels = ['linwcollapsed','linwcollapsed_cLOF','linwb','linwb_mrgLOF'], # kernels to consider 
        phenotype = lambda wc: phenotypes[ wc.pheno ],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/lrtsim_random/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        debug=False,
        random=True
    log:
        'logs/association/sclrt_kernels_missense_retest_random/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_missense_retest_top_hits.py'
        
        
rule assoc_missense_localcollapsing_retest_random_all:
    # run above rule for all phenotypes 
    input:
        expand(rules.assoc_missense_localcollapsing_retest_random.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])
        

rule assoc_spliceai_linw:
    # runs association tests for splice variants
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
        seak_install = rules.install_seak.output
    output:
        results_tsv = 'work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/results.tsv.gz',
        chi2param_tsv = 'work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/chi2mixparam.tsv.gz'
    params:
        phenotype = lambda wc: phenotypes[ wc.pheno ],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['splice_ai_min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/lrt_stats/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        debug = False,
        sclrt_nominal_significance_cutoff = 0.1
    log:
        'logs/association/sclrt_kernels_spliceai/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_spliceai.py'


rule assoc_spliceai_linw_all:
    input:
        expand(rules.assoc_spliceai_linw.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])

rule assoc_spliceai_linw_eval_top_hits:
    # calculates single-variant p-values and regression coefficients and other statistics for the most significant genes for splice variants
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
        seak_install = rules.install_seak.output
    output:
        touch('work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/top_hits/all.ok')
    params:
        kernels = ['linb','linwb','lin','linw','linb_mrgLOF','linwb_mrgLOF','lin_cLOF','linw_cLOF'],
        phenotype = lambda wc: phenotypes[ wc.pheno ],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['splice_ai_min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/top_hits/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        debug = False
    log:
        'logs/association/sclrt_kernels_spliceai_eval_top_hits/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_spliceai_eval_top_hits.py'

rule assoc_spliceai_linw_eval_top_hits_all:
    # runs rule above for all phenotypes - this rule will effectively run the entire pipeline for splice variants
    input:
         expand(rules.assoc_spliceai_linw_eval_top_hits.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])


rule assoc_spliceai_linw_retest_top_hits:
    # calculates single-variant p-values and regression coefficients and other statistics for the most significant genes for splice variants
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
        seak_install = rules.install_seak.output
    output:
        t = touch('work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/lrtsim/all.ok'),
        results_tsv = 'work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/lrtsim/lrt_retest.tsv.gz'
    params:
        kernels = ['linwb','linw','linwb_mrgLOF','linw_cLOF'],
        phenotype = lambda wc: phenotypes[wc.pheno],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['splice_ai_min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/lrtsim/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        debug = False,
        significance_cutoff = 5e-6,
        random=False
    log:
        'logs/association/sclrt_kernels_spliceai_retest_top_hits/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_spliceai_retest_top_hits.py'

rule assoc_spliceai_linw_retest_top_hits_all:
    # runs rule above for all phenotypes - this rule will effectively run the entire pipeline for splice variants
    input:
         expand(rules.assoc_spliceai_linw_retest_top_hits.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])

rule assoc_spliceai_linw_retest_random:
    # calculates single-variant p-values and regression coefficients and other statistics for the most significant genes for splice variants
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
        seak_install = rules.install_seak.output
    output:
        t = touch('work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/lrtsim_random/all.ok'),
        results_tsv = 'work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/lrtsim_random/lrt_retest.tsv.gz'
    params:
        kernels = ['linwb','linw','linwb_mrgLOF','linw_cLOF'],
        phenotype = lambda wc: phenotypes[wc.pheno],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['splice_ai_min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/lrtsim_random/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        debug = False,
        significance_cutoff = 5e-6,
        random = True
    log:
        'logs/association/sclrt_kernels_spliceai_retest_top_hits/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_spliceai_retest_top_hits.py'

rule assoc_spliceai_linw_retest_random_all:
    # runs rule above for all phenotypes - this rule will effectively run the entire pipeline for splice variants
    input:
         expand(rules.assoc_spliceai_linw_retest_random.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])

rule assoc_deepripe_single_localcollapsing:
    # run association tests with DeepRiPe variant effect predictions for single RBP
    input:
        bed = expand(rules.link_genotypes.output.bed, id = plinkfiles.getIds()),
        h5_rbp_plus = ancient(expand(rules.run_deepripe_vep.output.h5, id = plinkfiles.getIds(), strand=['plus'])),
        h5_rbp_minus = ancient(expand(rules.run_deepripe_vep.output.h5, id = plinkfiles.getIds(), strand=['minus'])),
        bed_rbp_plus = ancient(expand(rules.run_deepripe_vep.output.bed, id = plinkfiles.getIds(), strand=['plus'])),
        bed_rbp_minus = ancient(expand(rules.run_deepripe_vep.output.bed, id = plinkfiles.getIds(), strand=['minus'])),
        ensembl_vep_tsv = expand(rules.process_ensembl_vep_output_highimpact.output.tsv, id = plinkfiles.getIds(), allow_missing=True),
        mac_report = expand(rules.filter_variants.output.vid_tsv, id = plinkfiles.getIds(), allow_missing=True),
        regions_bed = rules.get_protein_coding_genes.output.pc_genes_bed,
        seak_install = rules.install_seak.output,
        covariates_tsv = config['covariates'],
        phenotypes_tsv = config['phenotypes']
    output:
        results_tsv = 'work/association/sclrt_kernels_deepripe/{filter_highconfidence}/{pheno}/{rbp}/results.tsv.gz',
        chi2param_tsv = 'work/association/sclrt_kernels_deepripe/{filter_highconfidence}/{pheno}/{rbp}/chi2mixparam.tsv.gz'
    params:
        phenotype = lambda wc: phenotypes[ wc.pheno ],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['deepripe_min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_deepripe/{filter_highconfidence}/{pheno}/{rbp}/lrt_stats/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno, rbp=wc.rbp),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        rbp_of_interest = lambda wc: [ wc.rbp ],
        debug = False
    log:
        'logs/association/sclrt_kernels_deepripe/{filter_highconfidence}_{rbp}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_deepripe_single.py'


rule assoc_deepripe_single_localcollapsing_all:
    # runs rule above for a list of RBPs of interest
    input:
        expand(rules.assoc_deepripe_single_localcollapsing.output, filter_highconfidence=['all'], rbp=['ELAVL1','HNRNPD','KHDRBS1_k562','MBNL1','QKI','QKI_hepg2','QKI_k562','TARDBP'], pheno=phenotypes.keys())

rule assoc_deepripe_single_localcollapsing_single:
    # runs rule above for a specific RBP
    input:
        expand(rules.assoc_deepripe_single_localcollapsing.output, filter_highconfidence=['all'], pheno=phenotypes.keys(), allow_missing=True)
    output:
        touch('work/association/sclrt_kernels_deepripe/{rbp}.ok')


rule assoc_deepripe_multiple_cholesky:
    # run association tests with DeepRiPe variant effect predictions for a group of RBPs
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
        phenotypes_tsv = config['phenotypes']
    output:
        results_tsv = 'work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/results.tsv.gz',
        chi2param_tsv = 'work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/chi2mixparam.tsv.gz'
    params:
        phenotype = lambda wc: phenotypes[ wc.pheno ],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['deepripe_min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/lrt_stats/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        rbp_of_interest = ['ELAVL1','HNRNPD','KHDRBS1_k562','MBNL1','QKI','QKI_hepg2','QKI_k562','TARDBP'],
        debug = False,
        sclrt_nominal_significance_cutoff = 0.1
    log:
        'logs/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}_{pheno}.log'
    resources:
        mem_mb=get_mem_mb(10000,1.5),
        time="8:00:00"
    threads:
        1
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_deepripe_multiple.py'


rule assoc_deepripe_multiple_cholesky_all:
    # runs rule above for all phenotypes
    input:
        expand(rules.assoc_deepripe_multiple_cholesky.output, filter_highconfidence=['all'], pheno=phenotypes.keys())

rule assoc_deepripe_multiple_cholesky_eval_top_hits:
    # calculates single-variant p-values and regression coefficients and other statistics for the most significant genes for DeepRiPe variants
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
        results_tsv = rules.assoc_deepripe_multiple_cholesky.output.results_tsv
    output:
        touch('work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/top_hits/all.ok')
    params:
        kernels = ['lincholesky', 'linwcholesky', 'lincholesky_notLOF', 'linwcholesky_notLOF'],
        phenotype = lambda wc: phenotypes[ wc.pheno ],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['deepripe_min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/top_hits/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        rbp_of_interest = rules.assoc_deepripe_multiple_cholesky.params.rbp_of_interest,
        debug = False
    resources:
        mem_mb=get_mem_mb(10000,1.5),
        time="01:00:00"
    threads:
        1
    log:
        'logs/association/sclrt_kernels_deepripe_multiple_eval_top_hits/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_deepripe_multiple_eval_top_hits.py'


rule assoc_deepripe_multiple_cholesky_eval_top_hits_all:
    # runs rule above for all phenotypes - this will effectively run the entire pipeline for DeepRiPe variants
    input:
        expand(rules.assoc_deepripe_multiple_cholesky_eval_top_hits.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])


rule assoc_deepripe_multiple_cholesky_retest_top_hits:
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
        results_tsv = rules.assoc_deepripe_multiple_cholesky.output.results_tsv
    output:
        t = touch('work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/lrtsim/all.ok'),
        results_tsv = 'work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/lrtsim/res_lrtsim.tsv.gz'
    params:
        significance_cutoff = 5e-6,
        kernels = ['linwcholesky', 'linwcholesky_notLOF'],
        phenotype = lambda wc: phenotypes[wc.pheno],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['deepripe_min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/lrtsim/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        rbp_of_interest = rules.assoc_deepripe_multiple_cholesky.params.rbp_of_interest,
        debug = False,
        random = False
    resources:
        mem_mb=get_mem_mb(10000,1.5),
        time="01:00:00"
    threads:
        1
    log:
        'logs/association/sclrt_kernels_deepripe_multiple_retest_top_hits/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_deepripe_multiple_retest_top_hits.py'


rule assoc_deepripe_multiple_cholesky_retest_top_hits_all:
    # runs rule above for all phenotypes - this will effectively run the entire pipeline for DeepRiPe variants
    input:
        expand(rules.assoc_deepripe_multiple_cholesky_retest_top_hits.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])
        

rule assoc_deepripe_multiple_cholesky_retest_random:
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
        results_tsv = rules.assoc_deepripe_multiple_cholesky.output.results_tsv
    output:
        t = touch('work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/lrtsim_random/all.ok'),
        results_tsv = 'work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/lrtsim_random/res_lrtsim.tsv.gz'
    params:
        significance_cutoff = 5e-6,
        kernels = ['linwcholesky', 'linwcholesky_notLOF'],
        phenotype = lambda wc: phenotypes[wc.pheno],
        covariate_column_names = config['covariate_column_names'],
        max_maf = config['maf_cutoff'],
        min_impact = config['deepripe_min_impact'],
        out_dir_stats = lambda wc: 'work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/lrtsim_random/'.format(filter_highconfidence=wc.filter_highconfidence, pheno=wc.pheno),
        ids = plinkfiles.getIds(),
        filter_highconfidence = lambda wc: {'all': False, 'highconf_only': True}[wc.filter_highconfidence],
        rbp_of_interest = rules.assoc_deepripe_multiple_cholesky.params.rbp_of_interest,
        debug = False,
        random = True
    log:
        'logs/association/sclrt_kernels_deepripe_multiple_retest_random/{filter_highconfidence}_{pheno}.log'
    conda:
        '../env/seak.yml'
    script:
        '../script/python/assoc_sclrt_kernels_deepripe_multiple_retest_top_hits.py'


rule assoc_deepripe_multiple_cholesky_retest_random_all:
    # runs rule above for all phenotypes - this will effectively run the entire pipeline for DeepRiPe variants
    input:
        expand(rules.assoc_deepripe_multiple_cholesky_retest_random.output, pheno=phenotypes.keys(), filter_highconfidence=['all'])

