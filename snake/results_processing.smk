


rule export_results_baseline:
    # export results and calculate genomic inflation
    input:
        rules.assoc_baseline_scoretest_all.input,
        rules.install_seak.output
    output:
        lambdaval='results/tables/lambdaval_baseline.tsv',
        results_tsv='results/tables/results_baseline.tsv.gz'
    conda:
        '../env/seak.yml'
    log:
        notebook='results/00_export_results_baseline.ipynb'
    notebook:
        '../notebooks/00_export_results_baseline.ipynb'
        
        
rule export_results_missense:
    # export results and calculate genomic inflation
    input:
        rules.assoc_missense_localcollapsing_all.input,
        rules.install_seak.output
    output:
        lambdaval='results/tables/lambdaval_missense.tsv',
        results_tsv='results/tables/results_missense.tsv.gz'
    conda:
        '../env/seak.yml'
    log:
        notebook='results/00_export_results_missense.ipynb'
    notebook:
        '../notebooks/00_export_results_missense.ipynb'
    

rule export_results_splice:
    # export results and calculate genomic inflation
    input:
        rules.assoc_spliceai_linw_all.input,
        rules.install_seak.output
    output:
        lambdaval='results/tables/lambdaval_splice.tsv',
        results_tsv='results/tables/results_splice.tsv.gz'
    conda:
        '../env/seak.yml'
    log:
        notebook='results/00_export_results_splice.ipynb'
    notebook:
        '../notebooks/00_export_results_splice.ipynb'


rule export_results_rbp:
    # export results and calculate genomic inflation
    input:
        rules.assoc_deepripe_multiple_cholesky_all.input,
        rules.install_seak.output
    output:
        lambdaval='results/tables/lambdaval_rbp.tsv',
        results_tsv='results/tables/results_rbp.tsv.gz'
    conda:
        '../env/seak.yml'
    log:
        notebook='results/00_export_results_rbp.ipynb'
    notebook:
        '../notebooks/00_export_results_rbp.ipynb'


rule lrt_pval_random_compare_missense:
    # compare gene-specific null distr. pvalues vs pooled (genome-wide) null distr. pvalues
    input:
        lrt_tsv=rules.assoc_missense_localcollapsing_retest_random.output['results_tsv'],
        chi2param=rules.assoc_missense_localcollapsing.output['chi2param_tsv'],
        install_seak_ok=rules.install_seak.output
    output:
        'work/association/sclrt_kernels_missense/{filter_highconfidence}/{pheno}/lrtsim_random/pval_comparison_lrt_pooled_vs_single.tsv.gz'
    conda:
        '../env/seak.yml'
    params:
        kern=['linwcollapsed','linwcollapsed_cLOF','linwb','linwb_mrgLOF']
    log:
        notebook='results/lrt_pooled_vs_single/{pheno}_{filter_highconfidence}_missense.ipynb'
    notebook:
        '../notebooks/01_lrt_pooled_vs_single_nulldistr_pval_compare.ipynb'
        
rule lrt_pval_random_compare_missense_all:
    input:
        expand(rules.lrt_pval_random_compare_missense.output, filter_highconfidence=['all'], pheno=phenotypes.keys())
    

rule lrt_pval_random_compare_splice:
    # compare gene-specific null distr. pvalues vs pooled (genome-wide) null distr. pvalues
    input:
        lrt_tsv=rules.assoc_spliceai_linw_retest_random.output['results_tsv'],
        chi2param=rules.assoc_spliceai_linw.output['chi2param_tsv'],
        install_seak_ok=rules.install_seak.output
    output:
        'work/association/sclrt_kernels_spliceai/{filter_highconfidence}/{pheno}/lrtsim_random/pval_comparison_lrt_pooled_vs_single.tsv.gz'
    conda:
        '../env/seak.yml'
    params:
        kern=['linw','linw_cLOF','linwb','linwb_mrgLOF']
    log:
        notebook='results/lrt_pooled_vs_single/{pheno}_{filter_highconfidence}_splice.ipynb'
    notebook:
        '../notebooks/01_lrt_pooled_vs_single_nulldistr_pval_compare.ipynb'

rule lrt_pval_random_compare_splice_all:
    input:
        expand(rules.lrt_pval_random_compare_splice.output, filter_highconfidence=['all'], pheno=phenotypes.keys())
        

rule lrt_pval_random_compare_rbp:
    # compare gene-specific null distr. pvalues vs pooled (genome-wide) null distr. pvalues
    input:
        lrt_tsv=rules.assoc_deepripe_multiple_cholesky_retest_random.output['results_tsv'],
        chi2param=rules.assoc_deepripe_multiple_cholesky.output['chi2param_tsv'],
        install_seak_ok=rules.install_seak.output
    output:
        'work/association/sclrt_kernels_deepripe_multiple/{filter_highconfidence}/{pheno}/lrtsim_random/pval_comparison_lrt_pooled_vs_single.tsv.gz'
    conda:
        '../env/seak.yml'
    params:
        kern=['linwcholesky']
    log:
        notebook='results/lrt_pooled_vs_single/{pheno}_{filter_highconfidence}_rbp.ipynb'
    notebook:
        '../notebooks/01_lrt_pooled_vs_single_nulldistr_pval_compare.ipynb'

rule lrt_pval_random_compare_rbp_all:
    input:
        expand(rules.lrt_pval_random_compare_rbp.output, filter_highconfidence=['all'], pheno=phenotypes.keys())
        
        
rule lrt_pval_random_compare_all:
    input:
        rules.lrt_pval_random_compare_missense_all.input,
        rules.lrt_pval_random_compare_splice_all.input,
        rules.lrt_pval_random_compare_rbp_all.input