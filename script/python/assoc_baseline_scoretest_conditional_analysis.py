


import logging
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)


import numpy as np
import pandas as pd
import gzip

from seak.scoretest import ScoretestNoK
from seak.data_loaders import intersect_ids, CovariatesLoaderCSV, VariantLoaderSnpReader
from pysnptools.snpreader import Bed
from util.association import BurdenLoaderHDF5
from util import Timer


'''
This script is the baseline that we later compare to.
It loops over genes and performs score-tests, using pre-computed indicator variables (0/1)
'''


# covariates loader
covariatesloader = CovariatesLoaderCSV(snakemake.params.phenotype,
                                       snakemake.input.covariates_tsv,
                                       snakemake.params.covariate_column_names,
                                       sep='\t',
                                       path_to_phenotypes=snakemake.input.phenotypes_tsv)

# set up burden loaders
bloader_lof = BurdenLoaderHDF5(snakemake.input.h5_lof, snakemake.input.iid_lof, snakemake.input.gid_lof)
bloader_missense = BurdenLoaderHDF5(snakemake.input.h5_missense, snakemake.input.iid_missense, snakemake.input.gid_missense)

# make sure individuals are in the same order
bloader_lof.update_individuals(covariatesloader.get_iids())
bloader_missense.update_individuals(covariatesloader.get_iids())

# gene names to iterate over
genes = np.union1d(bloader_lof.get_vids(), bloader_missense.get_vids())
if isinstance(genes, str):
    genes = [genes]
    
# keep only those that are significant!
# (1) load the results table
results = pd.read_csv(snakemake.input.results_tsv, sep='\t', na_values='.')

# the columns to consider
sigcols = snakemake.params.columns
if isinstance(sigcols, str):
    sigcols = [sigcols]
    
pvcols_score = sigcols
results = results[['gene', 'nCarrier_pLOF', 'nCarrier_missense'] + pvcols_score]
results['pheno'] = snakemake.params.phenotype
results = results.loc[results.gene.isin(genes)]
results['gene_name'] = results['gene'].str.split('_', expand=True)[1]
genes = results.gene.tolist()

# (2) load the list of variants to condition on
_conditional_list = pd.read_csv(snakemake.input.conditional_list, sep='\t', header=0, usecols=['gene_name','pheno'])
_conditional_list = _conditional_list[(_conditional_list.pheno == snakemake.wildcards['pheno']) | (_conditional_list.pheno == snakemake.params.phenotype) ]
_conditional_list = _conditional_list.drop_duplicates() 

if len(_conditional_list) == 0:
    logging.info('No genes pass significance threshold.'.format(snakemake.params.phenotype))
    with gzip.open(snakemake.output.results_tsv, 'wt') as outfile:
        outfile.write('# no significant hits for phenotype {}. \n'.format(snakemake.params.phenotype))
        sys.exit(0)
else:
    results =  results.merge(_conditional_list, how='right', on=['gene_name', 'pheno'])

# (3) keep only genes below threshold 
sig_genes = [results.gene[results[k] < snakemake.params.significance_cutoff].values for k in pvcols_score ]
sig_genes = np.unique(np.concatenate(sig_genes))

if len(sig_genes) == 0:
    logging.info('No genes pass significance threshold.')
    with gzip.open(snakemake.output.results_tsv, 'wt') as outfile:
        outfile.write('# no hits below specified threshold ({}) for phenotype {}. \n'.format(snakemake.params.significance_cutoff, snakemake.params.phenotype))
    sys.exit(0)
results = results.loc[results.gene.isin(sig_genes)].copy()
genes = results.gene.tolist()

logging.info('About to evaluate variants in {} genes'.format(len(results.gene.unique())))

# conditional analysis:
# these genotypes will be used for the conditional analysis
geno_cond = VariantLoaderSnpReader(Bed(snakemake.input.conditional_geno, count_A1=True, num_threads=1))
geno_cond.update_individuals(covariatesloader.get_iids())

# this file contains the mapping of associations to SNPs to condition on
conditional_list = pd.read_csv(snakemake.input.conditional_list, sep='\t', header=0)
conditional_list = conditional_list[(conditional_list.pheno == snakemake.params['phenotype']) | (conditional_list.pheno == snakemake.wildcards['pheno']) ].drop_duplicates()
conditional_list = conditional_list.loc[conditional_list.gene_name.isin(results.gene_name)]

geno_cond.update_variants(intersect_ids(conditional_list.snp_rsid, geno_cond.get_vids()))
logging.info('considering {} variants as covariates for conditional tests.'.format(len(geno_cond.get_vids())))

# set up the null model
Y, X = covariatesloader.get_one_hot_covariates_and_phenotype('NoK')
null_model = ScoretestNoK(Y, X)

logging.info('Phenotype: {}, Sample size: {}'.format(snakemake.params.phenotype, len(Y)))

def test_gene(gene):
    
    # conditional analysis
    # set up the function to load the variants to condition on
    
    def get_conditional(gene):
        
        # print(gene)
        # print(gene.split('_')[1])
        # print(conditional_list)
        
        cond_snps_vid = conditional_list.loc[ conditional_list.gene_name == gene.split('_')[1], 'snp_rsid' ].unique().tolist()
        
        temp_genotypes, temp_vids = geno_cond.genotypes_by_id(cond_snps_vid, return_pos=False)
        temp_genotypes -= np.nanmean(temp_genotypes, axis=0)
        
        cond_snps = np.ma.masked_invalid(temp_genotypes).filled(0.)

        return cond_snps, temp_vids
    
    pval_dict = {}
    pval_dict['gene'] = gene
    
    def call(GV, G2, name):
        # calls the test, appends results to pval_dict
        if GV is None:
            pval_dict['pv_' + name] = np.nan
            pval_dict['nCarrier_' + name] = 0
        else:
            pv = null_model.pv_alt_model(GV, G2)
            if pv < 0.:
                logging.warning('negative value encountered in p-value computation for gene {}, p-value: {}, using saddle instead.'.format(gene, pv))
                pv = null_model.pv_alt_model(GV, G2, method='saddle')
            pval_dict['pv_' + name] = pv
            pval_dict['nCarrier_' + name] = int(GV.sum())

            
    # Load the protein lof burdens
    try:
        Glof = bloader_lof.genotypes_by_id(gene).astype(float)
    except KeyError:
        Glof = None

    # Load the missense burdens
    try:
        Gmiss = bloader_missense.genotypes_by_id(gene).astype(float)
    except KeyError:
        Gmiss = None
        
    G2, cond_vid = get_conditional(gene)

    call(Glof, G2, 'pLOF')
    call(Gmiss, G2, 'missense')

    if (Glof is None) or (Gmiss is None):
        call(None, 'mrg')
    else:
        G = np.maximum(Glof, Gmiss)
        call(G, G2, 'mrg')
    
    pval_dict['cond_snps'] = ','.join(cond_vid)
    
    return pval_dict


timer = Timer()
results = pd.DataFrame.from_dict([test_gene(gene) for gene in genes]).set_index('gene')

t = timer.check()
logging.info('{} genes tested in {:.2f} minutes.'.format(len(results), t/60.))

results.to_csv(snakemake.output.results_tsv, sep='\t', index=True, na_rep='.')

