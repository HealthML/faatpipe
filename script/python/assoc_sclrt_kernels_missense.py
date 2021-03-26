


import os
# os.environ["OMP_NUM_THREADS"] = "16"

import logging
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)

import pandas as pd
import numpy as np

# seak imports
from seak.data_loaders import intersect_ids, EnsemblVEPLoader, VariantLoaderSnpReader, CovariatesLoaderCSV
from seak.kernels import LocalCollapsing
from seak.scoretest import ScoretestNoK
from seak.lrt import LRTnoK, pv_chi2mixture, fit_chi2mixture

from pysnptools.snpreader import Bed

from util.association import BurdenLoaderHDF5
from util import Timer


class GotNone(Exception):
    pass

# set up the covariatesloader

covariatesloader = CovariatesLoaderCSV(snakemake.params.phenotype,
                                       snakemake.input.covariates_tsv,
                                       snakemake.params.covariate_column_names,
                                       sep='\t',
                                       path_to_phenotypes=snakemake.input.phenotypes_tsv)

# initialize the null models
Y, X = covariatesloader.get_one_hot_covariates_and_phenotype('noK')

null_model_score = ScoretestNoK(Y, X)
null_model_lrt = LRTnoK(X, Y)


# set up function to filter variants:
def maf_filter(mac_report):

    # load the MAC report, keep only observed variants with MAF below threshold
    mac_report = pd.read_csv(mac_report, sep='\t', usecols=['SNP', 'MAF', 'Minor', 'alt_greater_ref'])
    
    if snakemake.params.filter_highconfidence:
        vids = mac_report.SNP[(mac_report.MAF < snakemake.params.max_maf) & (mac_report.Minor > 0) & ~(mac_report.alt_greater_ref.astype(bool)) & (mac_report.hiconf_reg.astype(bool))]
    else:
        vids = mac_report.SNP[(mac_report.MAF < snakemake.params.max_maf) & (mac_report.Minor > 0) & ~(mac_report.alt_greater_ref)] 

    # this has already been done in filter_variants.py
    # load the variant annotation, keep only variants in high-confidece regions
    # anno = pd.read_csv(anno_tsv, sep='\t', usecols=['Name', 'hiconf_reg'])
    # vids_highconf = anno.Name[anno.hiconf_reg.astype(bool).values]
    # vids = np.intersect1d(vids, vids_highconf)

    return mac_report.set_index('SNP').loc[vids]


# genotype path, vep-path:
assert len(snakemake.params.ids) == len(snakemake.input.bed), 'Error: length of chromosome IDs does not match length of genotype files'
geno_vep = zip(snakemake.params.ids, snakemake.input.bed, snakemake.input.vep_tsv, snakemake.input.mac_report, snakemake.input.h5_lof, snakemake.input.iid_lof, snakemake.input.gid_lof)

stats = []
simulations = []

i_gene = 0

# enter the chromosome loop:

timer = Timer()
for i, (chromosome, bed, vep_tsv, mac_report, h5_lof, iid_lof, gid_lof) in enumerate(geno_vep):
    
    
    # set up the ensembl vep loader for the chromosome
    ensemblvepdf = pd.read_csv(vep_tsv,
                               sep='\t',
                               usecols=['Uploaded_variation', 'Location', 'Gene', 'pos_standardized', 'impact'],
                               index_col='Uploaded_variation')

    # get set of variants for the chromosome:
    mac_report = maf_filter(mac_report)
    filter_vids = mac_report.index.values

    # filter by MAF
    keep = intersect_ids(filter_vids, ensemblvepdf.index.values)
    ensemblvepdf = ensemblvepdf.loc[keep]
    ensemblvepdf.reset_index(inplace=True)

    # filter by impact:
    ensemblvepdf = ensemblvepdf[ensemblvepdf.groupby(['Gene','pos_standardized'])['impact'].transform(np.max) >= snakemake.params.min_impact ]

    # initialize the loader
    eveploader = EnsemblVEPLoader(ensemblvepdf['Uploaded_variation'], ensemblvepdf['Location'], ensemblvepdf['Gene'], data = ensemblvepdf[['pos_standardized','impact']].values)

    # set up the regions to loop over for the chromosome
    regions = pd.read_csv(snakemake.input.regions_bed, sep='\t', header=None, usecols=[0,1,2,3], dtype={0:str, 1: np.int32, 2:np.int32, 3:str})
    regions.columns = ['chrom', 'start', 'end', 'name']

    # discard all genes for which we don't have annotations
    regions['gene'] = regions.name.str.split('_', expand=True)[0]
    regions.set_index('gene', inplace=True)
    genes = intersect_ids(np.unique(regions.index.values), np.unique(eveploader.pos_df.gene))
    regions = regions.loc[genes].reset_index()
    regions = regions.sort_values(['chrom','start','end'])[['chrom','start','end','name','gene']]

    # set up the variant loader (missense variants) for the chromosome
    plinkloader = VariantLoaderSnpReader(Bed(bed, count_A1=True, num_threads=4))
    plinkloader.update_variants(eveploader.get_vids())
    plinkloader.update_individuals(covariatesloader.get_iids())
    
    # set up the protein LOF burden loader
    bloader_lof = BurdenLoaderHDF5(h5_lof, iid_lof, gid_lof)
    bloader_lof.update_individuals(covariatesloader.get_iids())

    # set up local collapsing
    collapser = LocalCollapsing(distance_threshold=1.)

    # set up the missense genotype + vep loading function
    def get_missense(interval):

        try:
            V1 = eveploader.anno_by_interval(interval, gene=interval['name'].split('_')[0])
        except KeyError:
            raise GotNone

        if V1.index.empty:
            raise GotNone

        vids = V1.index.get_level_values('vid')
        V1 = V1.droplevel(['gene'])

        temp_genotypes, temp_vids = plinkloader.genotypes_by_id(vids, return_pos=False)

        temp_genotypes -= np.nanmean(temp_genotypes, axis=0)
        G1 = np.ma.masked_invalid(temp_genotypes).filled(0.)

        ncarrier = np.sum(G1 > 0.5, axis=0)
        cummac = mac_report.loc[vids].Minor
        
        # polyphen / sift impact
        weights = V1[1].values.astype(np.float64)

        # "standardized" positions -> codon start positions
        pos = V1[0].values.astype(np.int32)

        return G1, vids, weights, ncarrier, cummac, pos

    # set up the protein-LOF loading function

    def get_plof(interval):

        try:
            G2 = bloader_lof.genotypes_by_id(interval['name']).astype(np.float)
        except KeyError:
            G2 = None

        return G2

    # set up the test-function for a single gene
    def test_gene(interval, seed):

        pval_dict = {}
        sim_dict = {}

        pval_dict['gene'] = interval['name']
        sim_dict['gene'] = interval['name']

        def call_score(GV, name):
            pv = null_model_score.pv_alt_model(GV)
            if pv < 0.:
                logging.warning('negative value encountered in p-value computation for gene {}, p-value: {}, using saddle instead.'.format(interval['name'], pv))
                pv = null_model_score.pv_alt_model(GV, method='saddle')
            pval_dict['pv_score_' + name] = pv

        def call_lrt(GV, name):
            lik = null_model_lrt.altmodel(GV)
            sim = null_model_lrt.pv_sim(100, seed=seed)
            pval_dict['lrtstat_' + name] = lik['stat']
            pval_dict['alteqnull_' + name] = float(lik['alteqnull'])
            sim_dict[name] = sim['res']

        # load missense variants
        G1, vids, weights, ncarrier, cummac, pos = get_missense(interval)

        # do a score burden test (max weighted), this is different than the baseline!
        G1_burden = np.max(np.where(G1 > 0.5, np.sqrt(weights), 0.), axis=1, keepdims=True)
        call_score(G1_burden, 'linwb')

        # perform local collapsing with weights
        if G1.shape[1] > 1:
            G1, clusters = collapser.collapse(G1, pos, np.sqrt(weights)) # will lead to a crash when there are negative weights...
        else:
            G1 = G1.dot(np.diag(np.sqrt(weights), k=0))
            clusters = [0]

        # do a score test (local collapsing)
        call_score(G1, 'linwcollapsed')

        # if gene is nominally significant:
        # if (pval_dict['pv_score_linwb'] < 0.01) | (pval_dict['pv_score_linwcollapsed'] < 0.01):
        if (pval_dict['pv_score_linwb'] < snakemake.params.sclrt_nominal_significance_cutoff) | (pval_dict['pv_score_linwcollapsed'] < snakemake.params.sclrt_nominal_significance_cutoff):

            # do lrt tests 
            call_lrt(G1_burden, 'linwb')
            call_lrt(G1, 'linwcollapsed')

            # load plof burden
            G2 = get_plof(interval)

            if G2 is not None:

                # merged (single variable)
                G1_burden_mrg = np.maximum(G2, G1_burden)
                call_score(G1_burden_mrg, 'linwb_mrgLOF')
                call_lrt(G1_burden_mrg, 'linwb_mrgLOF')

                # concatenated ( >= 2 variables)
                G1 = np.concatenate([G1, G2], axis=1)
                call_score(G1, 'linwcollapsed_cLOF')
                call_lrt(G1, 'linwcollapsed_cLOF')

        pval_dict['nCarrier'] = ncarrier.sum()
        pval_dict['cumMAC'] = cummac.sum()
        pval_dict['n_snp'] = len(vids)
        pval_dict['n_cluster'] = len(set(clusters))

        # sanity check
        assert pval_dict['cumMAC'] >= pval_dict['nCarrier'], 'Error: something is broken.'
        
        return pval_dict, sim_dict

    logging.info('loaders for chromosome {} initialized in {:.1f} seconds.'.format(chromosome, timer.check()))
    # run tests for all genes on the chromosome
    for _, region in regions.iterrows():

        try:
            gene_stats, sims = test_gene(region, i_gene)
        except GotNone:
            continue

        stats.append(gene_stats)
        simulations.append(sims)
        i_gene += 1
        if (i_gene % 100) == 0:
            logging.info('tested {} genes...'.format(i_gene))
        # print(i_gene)

    logging.info('all tests for chromosome {} performed in {:.2f} minutes.'.format(chromosome, timer.check()/60.))
    logging.info('genes tested so far: {}'.format(i_gene + 1))


# when all chromosomes are done:

# generate results table:
stats = pd.DataFrame.from_dict(stats).set_index('gene')

# these are the kernels 
kern = ['linwb', 'linwb_mrgLOF', 'linwcollapsed', 'linwcollapsed_cLOF']

# concatenate and save gene-specific simulations
# ...complicated nested dict comprehension that ensures we drop empty values etc.
simulations_ = { k: {gene: simval for gene, simval in ((s['gene'], s[k]) for s in simulations if k in s) if len(simval) > 0} for k in kern}

# print(simulations_)

if not os.path.isdir(snakemake.params.out_dir_stats):
    os.makedirs(snakemake.params.out_dir_stats)

for k in kern:
    pd.DataFrame.from_dict(simulations_[k], orient='index').to_pickle(snakemake.params.out_dir_stats + '{}.pkl.gz'.format(k))

# calculate chi2 mixture parameters for each kernel
params = {k : fit_chi2mixture(np.concatenate(list(simulations_[k].values())), 0.1) for k in simulations_.keys()}

pvals = np.empty((stats.shape[0], len(kern)))
pvals[:, :] = np.nan
pvals = pd.DataFrame(pvals, index=stats.index.values, columns=[ 'pv_lrt_' + k for k in kern ])

for k in simulations_.keys():
    mask = ~np.isnan(stats['lrtstat_' + k].values)
    pvals.loc[mask, 'pv_lrt_' + k] = pv_chi2mixture(stats.loc[mask, 'lrtstat_' + k].values, params[k]['scale'], params[k]['dof'], params[k]['mixture'], alteqnull=stats.loc[mask, 'alteqnull_' + k].values.astype(bool))

pvals = stats.join(pvals, how='left')
pvals['pheno'] = snakemake.params.phenotype

pvals.to_csv(snakemake.output.results_tsv, sep='\t', index=True)

params = pd.DataFrame(params)
params.to_csv(snakemake.output.chi2param_tsv, sep='\t', index=True)


