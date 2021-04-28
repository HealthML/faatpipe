
# These rules implement variant effect prediction using DeepRiPe
# the {id} wildcard refers to the different chromosomes

rule split_vcf_strand:
    # DeepRiPe produces strand-specific predictions
    # We split the VCF into variants overlapping genes on the plus/minus strand
    # produces tabix-indexed vcf files corresponding to variants on the separate strands that overlap protein-coding genes
    input:
        pc_genes_bed_plus = rules.get_protein_coding_genes.output.pc_genes_bed_plus,
        pc_genes_bed_minus = rules.get_protein_coding_genes.output.pc_genes_bed_minus,
        vcf = rules.bim_to_vcf.output.vcf
    output:
        vcf_plus = 'work/variant_effect_prediction/vcf/by_strand/{id}_plus.recode.vcf.gz',
        vcf_minus = 'work/variant_effect_prediction/vcf/by_strand/{id}_minus.recode.vcf.gz'
    params:
        out_prefix_plus = lambda wc, output: output.vcf_plus.replace('.recode.vcf.gz',''),
        out_prefix_minus = lambda wc, output: output.vcf_minus.replace('.recode.vcf.gz',''),
        vcf_plus_tmp = lambda wc, output: output.vcf_plus.replace('.gz',''),
        vcf_minus_tmp = lambda wc, output: output.vcf_minus.replace('.gz','')
    conda:
        '../env/genomics.yml'
    log:
        'logs/deepripe/split_vcf_strand_{id}.log'
    shell:
        "("
        "vcftools --bed {input.pc_genes_bed_plus} --gzvcf {input.vcf} --out {params.out_prefix_plus} --recode; "
        "vcftools --bed {input.pc_genes_bed_minus} --gzvcf {input.vcf} --out {params.out_prefix_minus} --recode; "
        "bgzip {params.vcf_plus_tmp} && bgzip {params.vcf_minus_tmp}; "
        "tabix {output.vcf_plus} && tabix {output.vcf_minus} "
        ") &> {log}"


rule merge_vcf_strand:
    # Merge the strand-specific VCF files
    # this will break if the plink-file identifiers (which are assumed to be chromosomes) don't follow chromosomal order!
    # -> currently not needed / dead end
    input:
        vcf = expand('work/variant_effect_prediction/vcf/by_strand/{id}_{{strand}}.recode.vcf.gz', id=sorted(plinkfiles.getIds()))
    output:
        vcf_mrg = 'work/variant_effect_prediction/vcf/by_strand/mrg/{strand}.vcf.gz',
        vcf_mrg_index = 'work/variant_effect_prediction/vcf/by_strand/mrg/{strand}.vcf.gz.tbi'
    conda:
        '../env/genomics.yml'
    params:
        vcf = lambda wc, input: ' '.join(input.vcf)
    log:
        'logs/deepripe/mrg_vcf_{strand}.log'
    shell:
        '('
        'bcftools concat {params.vcf} -Oz -o {output.vcf_mrg}; '
        'tabix {output.vcf_mrg}; '
        ') &> {log}'
        
        
rule run_deepripe_vep:
    # runs the DeepRiPe variant effect prediction for all chromosomes + strands independently
    # this rule requires a GPU
    input:
        vcf = 'work/variant_effect_prediction/vcf/by_strand/{id}_{strand}.recode.vcf.gz',
        ref_fa = config['reference_fa']
    output:
        h5 = 'work/variant_effect_prediction/deepripe/{id}/{strand}_scores.h5',
        bed = 'work/variant_effect_prediction/deepripe/{id}/{strand}.bed.gz',
        variants_txt = 'work/variant_effect_prediction/deepripe/{id}/{strand}.ids.txt' 
    params:
        out_prefix = 'work/variant_effect_prediction/deepripe/{id}/{strand}'
    log:
        'logs/deepripe/run_deepripe_vep_{id}_{strand}.log'
    conda:
        '../env/deepripe.yml'
    resources:
        gpu=0
    shell:
        "("
        "python script/python/run_vep_deepripe.py "
        "-gpu {resources.gpu} "
        "-strand {wildcards.strand} "
        "-vcf {input.vcf} "
        "-ref {input.ref_fa} "
        "-out {params.out_prefix} "
        ") &> {log}"
        

rule run_deepripe_vep_all:
    # rule to run all DeepRiPe predictions
    input:
        expand(rules.run_deepripe_vep.output, id=plinkfiles.getIds(), strand=['plus', 'minus'])
        
