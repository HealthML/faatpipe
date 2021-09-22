
import os
import sys
# os.environ["OMP_NUM_THREADS"] = "16"

import logging
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)

import pandas as pd
import numpy as np
import pickle

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

def get_regions():
    # load the results, keep those below a certain p-value
    results = pd.read_csv(snakemake.input.results_tsv, sep='\t')

    kern = snakemake.params.kernels
    if isinstance(kern, str):
        kern = [kern]

    pvcols_score = ['pv_score_' + k for k in kern ]
    pvcols_lrt = ['pv_lrt_' + k for k in kern]
    statcols = ['lrtstat_' + k for k in kern]
    results = results[['gene', 'n_snp', 'cumMAC', 'nCarrier'] + statcols + pvcols_score + pvcols_lrt]

    # get genes below threshold
    genes = [results.gene[results[k] < 1e-7].values for k in pvcols_score + pvcols_lrt ]
    genes = np.unique(np.concatenate(genes))

    if len(genes) == 0:
        return None

    # set up the regions to loop over for the chromosome
    regions = pd.read_csv(snakemake.input.regions_bed, sep='\t', header=None, usecols=[0 ,1 ,2 ,3, 5], dtype={0 :str, 1: np.int32, 2 :np.int32, 3 :str, 5:str})
    regions.columns = ['chrom', 'start', 'end', 'name', 'strand']
    regions['strand'] = regions.strand.map({'+': 'plus', '-': 'minus'})
    regions = regions.set_index('name').loc[genes]

    regions = regions.join(results.set_index('gene'), how='left').reset_index()
    return regions


# genotype path, vep-path:
assert len(snakemake.params.ids) == len(snakemake.input.bed), 'Error: length of chromosome IDs does not match length of genotype files'
geno_vep = zip(snakemake.params.ids, snakemake.input.bed, snakemake.input.vep_tsv, snakemake.input.mac_report, snakemake.input.h5_lof, snakemake.input.iid_lof, snakemake.input.gid_lof)


regions_all = get_regions()
if regions_all is None:
    logging.info('No genes pass significance threshold, exiting.')
    sys.exit(0)

logging.info('About to evaluate variants in {} genes'.format(len(regions_all)))

i_gene = 0
# where we store the results
stats = []
# enter the chromosome loop:
timer = Timer()
for i, (chromosome, bed, vep_tsv, mac_report, h5_lof, iid_lof, gid_lof) in enumerate(geno_vep):

    if chromosome.replace('chr','') not in regions_all.chrom.unique():
        continue

    # set up the ensembl vep loader for the chromosome
    ensemblvepdf = pd.read_csv(vep_tsv,
                               sep='\t',
                               usecols=['Uploaded_variation', 'Location', 'Gene', 'pos_standardized', 'impact', 'ref', 'alt', 'cosine_similarity'],
                               index_col='Uploaded_variation')

    # get set of variants for the chromosome:
    mac_report = maf_filter(mac_report)
    filter_vids = mac_report.index.values

    # filter by MAF
    keep = intersect_ids(filter_vids, ensemblvepdf.index.values)
    ensemblvepdf = ensemblvepdf.loc[keep]
    ensemblvepdf.reset_index(inplace=True)

    # filter by impact:
    ensemblvepdf = ensemblvepdf[ensemblvepdf.groupby(['Gene' ,'pos_standardized'])['impact'].transform(np.max) >= snakemake.params.min_impact ]

    # initialize the loader
    eveploader = EnsemblVEPLoader(ensemblvepdf['Uploaded_variation'], ensemblvepdf['Location'], ensemblvepdf['Gene'], data = ensemblvepdf[['pos_standardized' ,'impact', 'ref', 'alt', 'cosine_similarity']].values)

    # set up the regions to loop over for the chromosome
    regions = regions_all.copy()

    # discard all genes for which we don't have annotations
    regions['gene'] = regions.name.str.split('_', expand=True)[0]
    regions.set_index('gene', inplace=True)
    genes = intersect_ids(np.unique(regions.index.values), np.unique(eveploader.pos_df.gene))
    regions = regions.loc[genes].reset_index()
    regions = regions.sort_values(['chrom' ,'start' ,'end'])

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

        ncarrier = np.nansum(np.nansum(temp_genotypes, axis=1) >= 1)

        temp_genotypes -= np.nanmean(temp_genotypes, axis=0)
        G1 = np.ma.masked_invalid(temp_genotypes).filled(0.)

        cummac = mac_report.loc[vids].Minor

        # polyphen / sift impact
        weights = V1[1].values.astype(np.float64)

        # "standardized" positions -> codon start positions
        pos = V1[0].values.astype(np.int32)

        ref = V1[2].values.astype(str)
        alt = V1[3].values.astype(str)
        cosine_similarity = V1[4].values.astype(np.float64)

        return G1, vids, weights, ncarrier, cummac, pos, ref, alt, cosine_similarity

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
        pval_dict['gene'] = interval['name']
        called = []

        def pv_score(GV):
            pv = null_model_score.pv_alt_model(GV)
            if pv < 0.:
                pv = null_model_score.pv_alt_model(GV, method='saddle')
            return pv

        def call_score(GV, name, vids=None):
            if name not in pval_dict:
                pval_dict[name] = {}
                called.append(name)
            pval_dict[name] = {}
            # single-marker p-values
            pval_dict[name]['pv_score'] = np.array([pv_score(GV[:,i,np.newaxis]) for i in range(GV.shape[1])])

            # single-marker coefficients 
            beta = [ null_model_score.coef(GV[:,i,np.newaxis]) for i in range(GV.shape[1]) ]
            pval_dict[name]['beta'] = np.array([x['beta'][0,0] for x in beta]) 
            pval_dict[name]['betaSd'] = np.array([np.sqrt(x['var_beta'][0,0]) for x in beta])
            if vids is not None:
                pval_dict[name]['vid'] = vids


        def call_lrt(GV, name, vids=None):
            if name not in pval_dict:
                pval_dict[name] = {}
                called.append(name)
            # get gene parameters, test statistics and and single-marker regression weights
            lik = null_model_lrt.altmodel(GV)
            pval_dict[name]['nLL'] = lik['nLL']
            pval_dict[name]['sigma2'] = lik['sigma2']
            pval_dict[name]['lrtstat'] = lik['stat']
            pval_dict[name]['h2'] = lik['h2']
            logdelta = null_model_lrt.model1.find_log_delta(GV.shape[1])
            pval_dict[name]['log_delta'] = logdelta['log_delta']
            pval_dict[name]['coef_random'] = null_model_lrt.model1.getPosteriorWeights(logdelta['beta'], logdelta=logdelta['log_delta'])
            if vids is not None:
                pval_dict[name]['vid'] = vids


        # load missense variants
        G1, vids, weights, ncarrier, cummac, pos, ref, alt, cosine_similarity = get_missense(interval)

        # these are common to all kernels
        pval_dict['vid'] = vids
        pval_dict['weights'] = weights
        pval_dict['MAC'] = cummac
        pval_dict['nCarrier'] = ncarrier
        pval_dict['pos'] = pos
        pval_dict['ref'] = ref
        pval_dict['alt'] = alt
        pval_dict['cosine_similarity'] = cosine_similarity

        # single-variant p-values:
        call_score(G1, 'variant_pvals') # score-pvalues and coeficients estimated independently
        call_lrt(G1.dot(np.diag(np.sqrt(weights), k=0)), 'variant_pvals') # coeficients after weight adjustment (estimated *jointly*)

        # sanity checks
        assert len(vids) == interval['n_snp'], 'Error: number of variants does not match! expected: {}  got: {}'.format(interval['n_snp'], len(vids))
        assert cummac.sum() == interval['cumMAC'], 'Error: cumMAC does not match! expeced: {}, got: {}'.format(interval['cumMAC'], cummac.sum())

        # do a score burden test (max weighted), this is different than the baseline!
        G1_burden = np.max(np.where(G1 > 0.5, np.sqrt(weights), 0.), axis=1, keepdims=True)

        # burden
        call_score(G1_burden, 'linwb')
        call_lrt(G1_burden, 'linwb')

        # perform local collapsing with weights
        if G1.shape[1] > 1:
            G1, clusters = collapser.collapse(G1, pos, np.sqrt(weights)) # will lead to a crash when there are negative weights...
        else:
            G1 = G1.dot(np.diag(np.sqrt(weights), k=0))
            clusters = np.array([0])

        pval_dict['cluster_id'] = clusters

        # get the single-cluster p-values + stats
        call_score(G1, 'linwcollapsed')
        call_lrt(G1, 'linwcollapsed')
        pval_dict['linwcollapsed']['cluster_id'] = sorted(set(clusters))

        # load plof burden
        G2 = get_plof(interval)

        if G2 is not None:

            call_score(G2, 'LOF')
            call_lrt(G2, 'LOF')

            # merged (single variable)
            G1_burden_mrg = np.maximum(G2, G1_burden)
            call_score(G1_burden_mrg, 'linwb_mrgLOF')
            call_lrt(G1_burden_mrg, 'linwb_mrgLOF')

            call_score(np.concatenate([G1, G2], axis=1), 'linwcollapsed_cLOF')
            call_lrt(np.concatenate([G1, G2], axis=1), 'linwcollapsed_cLOF')
            pval_dict['linwcollapsed_cLOF']['cluster_id'] = sorted(set(clusters)) + [-1] # -1 indicates the LOF cluster

        return pval_dict, called

    logging.info('loaders for chromosome {} initialized in {:.1f} seconds.'.format(chromosome, timer.check()))
    # run tests for all genes on the chromosome
    for _, region in regions.iterrows():

        try:
            gene_stats, called = test_gene(region, i_gene)
        except GotNone:
            continue


        # build the single-variant datafame
        single_var_columns = ['gene','vid','weights','MAC','nCarrier','pos','ref','alt','cosine_similarity','cluster_id']
        sv_df = pd.DataFrame.from_dict({k: gene_stats[k] for k in single_var_columns})

        sv_df['pv_score'] = gene_stats['variant_pvals']['pv_score'] # single-variant p-values estimated independently
        sv_df['coef_random'] = gene_stats['variant_pvals']['coef_random']  # single-variant coefficients estimated jointly after weighting
        sv_df['beta'] = gene_stats['variant_pvals']['beta'] # single-variant coeffcients estimated independently *without* weighting
        sv_df['betaSd'] = gene_stats['variant_pvals']['betaSd'] # standard errors for the single-variant coefficients estimated independently *without* weighting
        sv_df['pheno'] = snakemake.params.phenotype

        out_dir = os.path.join(snakemake.params.out_dir_stats, region['name'])
        os.makedirs(out_dir, exist_ok=True)

        sv_df.to_csv(out_dir + '/variants.tsv.gz', sep='\t', index=False)


        for k in called:

            if k == 'variant_pvals':
                continue

            results_dict = gene_stats[k]
                
            df_cols = ['cluster_id','pv_score', 'coef_random', 'beta', 'betaSd', 'vid'] # parts of the dict that have lenght > 1
            
            df = pd.DataFrame.from_dict(data={k: results_dict[k] for k in df_cols if k in results_dict})
            df['gene'] = gene_stats['gene']
            df['pheno'] = snakemake.params.phenotype
            df.to_csv(out_dir + '/{}.tsv.gz'.format(k), sep='\t', index=False)
            
            #  other cols ['nLL', 'sigma2', 'lrtstat', 'h2', 'log_delta']
            other_cols = {k: v for k, v in results_dict.items() if k not in df_cols }
            other_cols['gene'] = gene_stats['gene']
            other_cols['pheno'] = snakemake.params.phenotype
            
            pickle.dump(other_cols, open(out_dir + '/{}_stats.pkl'.format(k), 'wb'))
            
                            

        i_gene += 1
        logging.info('tested {} genes...'.format(i_gene))
        
    timer.reset()

