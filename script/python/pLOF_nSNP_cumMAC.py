import pandas as pd
import numpy as np

pheno = snakemake.wildcards.pheno

filter_hiconf = {'all': False, 'highconf_only': True}[snakemake.wildcards.filter_highconfidence]

results = []

genes = pd.read_csv(snakemake.input.pc_genes, sep='\t', header=None, dtype=str)

genes['ensembl_id'] = genes[3].str.split('_',expand=True)[0]

for chrom in range(1,23):
    
    
    hc = pd.read_csv('work/basic_annotation/chr{}.tsv.gz'.format(chrom), sep='\t', usecols=['Name', 'hiconf_reg'])
    mr = pd.read_csv('work/mac_report/{}/chr{}.tsv.gz'.format(pheno, chrom), sep='\t', header=0)
    
    anno = pd.read_csv('work/variant_effect_prediction/ensembl_vep/processed/high_impact/chr{}.tsv.gz'.format(chrom), sep='\t')
    
    anno = anno.merge(mr, how='inner', left_on='Uploaded_variation', right_on='SNP', validate='many_to_one')
    anno = anno.merge(hc, how='left', left_on='Uploaded_variation', right_on='Name', validate='many_to_one')
    
    anno = anno[~anno.alt_greater_ref.astype(bool)]
    anno = anno[anno.MAF <= snakemake.config['maf_cutoff']]
    
    if filter_hiconf:
        anno = anno[anno.hiconf_reg.astype(bool)]
    
    groups = anno.groupby('Gene')['Uploaded_variation'].groups
    
    def get_cummac(gene):
        try:
            return anno.loc[groups[gene]].Minor.sum()
        except KeyError:
            return 0
        
    def get_n_snp(gene):
        try:
            return len(groups[gene])
        except KeyError:
            return 0
    
    genes_chr = genes[genes[0] == str(chrom)].copy()
    
    genes_chr['n_snp_pLOF'] = genes_chr.ensembl_id.transform(get_n_snp)
    genes_chr['cumMAC_pLOF'] = genes_chr.ensembl_id.transform(get_cummac)
    
    results.append(genes_chr)
    
    #print(chrom)

genes = pd.concat(results)[['ensembl_id', 'n_snp_pLOF', 'cumMAC_pLOF']]

genes.to_csv(snakemake.output.tsv, sep='\t', index=False)