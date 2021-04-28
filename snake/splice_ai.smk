
# These rules implement overlapping UK Biobank variants with SpliceAI variant effect predictions 
# the {id} wildcard refers to the different chromosomes

rule splice_ai_tabix_vcf:
    # index the vcf file with tabix
    input:
        vcf = config['splice_ai_scores_snv_vcf']
    output:
        vcf_tbi = config['splice_ai_scores_snv_vcf']+'.tbi'
    log:
        'logs/splice_ai_tabix_vcf.log'
    conda:
        '../env/genomics.yml'
    shell:
        '('
        'tabix {input.vcf}'
        ') &> {log}'


rule splice_ai_vcf_to_tsv:
    # converts the Splice AI VCF to a TSV file, keeps only those positions with effects >= 0.1
    # this GREATLY reduces the file size
    input:
        vcf = config['splice_ai_scores_snv_vcf'],
        vcf_tbi = rules.splice_ai_tabix_vcf.output.vcf_tbi
    output:
        tsv = 'data/splice_ai/splice_ai_scores_snv_{id}.tsv.gz'
    log:
        'logs/splice_ai/splice_ai_vcf_to_tsv_{id}.log'
    conda:
        '../env/genomics.yml'
    shell:
        "("
        "bash script/bash/splice_ai_vcf_to_tsv.sh {input.vcf} {wildcards.id} | gzip > {output.tsv};"
        ") &> {log}"


rule splice_ai_vcf_to_tsv_all:
    # runs the rule above for all chromosomes
    input:
        expand(rules.splice_ai_vcf_to_tsv.output, id=plinkfiles.getIds())

        
rule splice_ai_filter_and_overlap_with_genotypes:
    # overlaps the genotype calls in a bim file (i.e. chromosome) with the filtered spliceAI predictions
    input:
        tsv = rules.splice_ai_vcf_to_tsv.output.tsv,
        bim = lambda wc : plinkfiles.get('bim', wc.id)
    output:
        tsv = 'work/variant_effect_prediction/splice_ai/splice_ai_snv_{id}.tsv.gz'
    log:
        'logs/splice_ai/splice_ai_filter_and_overlap_with_genotypes_{id}.log'
    conda:
        '../env/genomics.yml'
    threads:
        1
    script:
        '../script/python/splice_ai_filter_and_overlap_with_genotypes.py'


rule splice_ai_filter_and_overlap_with_genotypes_all:
    # runs the rule above for all chromosomes
    input:
        expand(rules.splice_ai_filter_and_overlap_with_genotypes.output, id = plinkfiles.getIds())
    threads:
        1


rule splice_ai_merge_tsv:
    # simply concatenates the outputs of the rule above for all bim files
    # currently not needed / dead end
    input:
        tsv = expand(rules.splice_ai_filter_and_overlap_with_genotypes.output.tsv, id = sorted(plinkfiles.getIds()))
    output:
        tsv = protected('work/splice_ai/mrg/splice_ai_snv.tsv.gz')
    params:
        input_tsv = lambda wc, input : ' '.join(input.tsv)
    log:
        'logs/splice_ai/splice_ai_merge_tsv.log'
    run:
        for i, file in enumerate(input.tsv):
            if i == 0:
                shell('(cp {file} {{output.tsv}}) &> {{log}}'.format(file=file))
            else:
                # if we didn't have the header line we wouldn't have to do it this way...
                shell("(zcat {file} | awk 'NR>1' | gzip >> {{output.tsv}}) &> {{log}}".format(file=file))
        

