


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
        lambdaval='results/tables/lambdaval_rbp.tsv',
        results_tsv='results/tables/results_rbp.tsv.gz'
    conda:
        '../env/seak.yml'
    log:
        notebook='results/00_export_results_splice.ipynb'
    notebook:
        '../notebooks/00_export_results_splice.ipynb'


rule export_results_rbp:
    # export results and calculate genomic inflation
    input:
        rules.assoc_deepripe_single_localcollapsing_all.input,
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




