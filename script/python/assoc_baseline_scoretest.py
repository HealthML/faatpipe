


import logging
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)


import numpy as np
import pandas as pd

from seak.scoretest import ScoretestNoK
from seak.data_loaders import intersect_ids, CovariatesLoaderCSV
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

# set up the null model
Y, X = covariatesloader.get_one_hot_covariates_and_phenotype('NoK')
null_model = ScoretestNoK(Y, X)

logging.info('Phenotype: {}, Sample size: {}'.format(snakemake.params.phenotype, len(Y)))

def test_gene(gene):

    pval_dict = {}

    pval_dict['gene'] = gene
    
    def call(GV, name):
        # calls the test, appends results to pval_dict
        if GV is None:
            pval_dict['pv_' + name] = np.nan
            pval_dict['nCarrier_' + name] = 0
            pval_dict['beta_' + name ] = np.nan
            pval_dict['betaSd_' + name ] = np.nan
        else:
            pv = null_model.pv_alt_model(GV)
            if pv < 0.:
                logging.warning('negative value encountered in p-value computation for gene {}, p-value: {}, using saddle instead.'.format(gene, pv))
                pv = null_model.pv_alt_model(GV, method='saddle')
            pval_dict['pv_' + name] = pv
            pval_dict['nCarrier_' + name] = int(GV.sum())
            if pv < 1e-5:
                # if gene is quite significant get the regression coefficient + SE
                beta = null_model.coef(GV)
                pval_dict['beta_' + name ] = beta['beta'][0,0]
                pval_dict['betaSd_' + name ] = np.sqrt(beta['var_beta'][0,0])

    # Load the protein lof burdens
    try:
        Glof = bloader_lof.genotypes_by_id(gene).astype(np.float)
    except KeyError:
        Glof = None

    # Load the missense burdens
    try:
        Gmiss = bloader_missense.genotypes_by_id(gene).astype(np.float)
    except KeyError:
        Gmiss = None

    call(Glof, 'pLOF')
    call(Gmiss, 'missense')

    if (Glof is None) or (Gmiss is None):
        call(None, 'mrg')
    else:
        G = np.maximum(Glof, Gmiss)
        call(G, 'mrg')
            
    return pval_dict


timer = Timer()
results = pd.DataFrame.from_dict([test_gene(gene) for gene in genes]).set_index('gene')

t = timer.check()
logging.info('{} genes tested in {:.2f} minutes.'.format(len(results), t/60.))

results.to_csv(snakemake.output.results_tsv, sep='\t', index=True, na_rep='.')

