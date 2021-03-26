import os
# os.environ["OMP_NUM_THREADS"] = "16"

import logging

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)

import pandas as pd
import numpy as np
import h5py
import gzip
import pickle
import sys

from numpy.linalg import cholesky, LinAlgError
from sklearn.metrics.pairwise import rbf_kernel, cosine_similarity

# seak imports
from seak.data_loaders import intersect_ids, Hdf5Loader, VariantLoaderSnpReader, CovariatesLoaderCSV
from seak.kernels import LocalCollapsing
from seak.scoretest import ScoretestNoK
from seak.lrt import LRTnoK, pv_chi2mixture, fit_chi2mixture

from pysnptools.snpreader import Bed

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
        # note: We don't filter out the variants for which alt/ref are "flipped" for the RBPs
        vids = mac_report.SNP[(mac_report.MAF < snakemake.params.max_maf) & (mac_report.Minor > 0) & (mac_report.hiconf_reg.astype(bool))]
    else:
        vids = mac_report.SNP[(mac_report.MAF < snakemake.params.max_maf) & (mac_report.Minor > 0)]

    return mac_report.set_index('SNP').loc[vids]


def vep_filter(h5_rbp, bed_rbp):
    # returns variant ids for variants that pass filtering threshold and mask for the hdf5loader

    with h5py.File(h5_rbp, 'r') as f:

        rbp_of_interest = snakemake.params.rbp_of_interest
        try:
            # in newer version of h5py the strings are not decoded -> decode them
            labels = [l.decode() for l in f['labels'][:]]
        except AttributeError:
            # this used to work in older versions of  h5py:
            labels = list(f['labels'][:])

        if isinstance(rbp_of_interest, str):
            rbp_of_interest = [rbp_of_interest]

        assert all([rbp in labels for rbp in rbp_of_interest]), 'Error: not all of {} are in labels: {}'.format(rbp_of_interest, labels)

        keep = list(i for i, l in enumerate(labels) if l in rbp_of_interest)
        diffscores = f['diffscore'][:, keep]

    rbp_mask = np.array([i in rbp_of_interest for i in labels])

    if np.ndim(diffscores) == 1:
        diffscores = diffscores[:, np.newaxis]

    vid_mask = np.max(np.abs(diffscores), axis=1) >= snakemake.params.min_impact

    vid_df = pd.read_csv(bed_rbp, sep='\t', header=None, usecols=[3], names=['vid'], dtype={'vid': str})
    vids = vid_df.vid.str.split(pat='_[ACGT]+>[ACGT]+$', n=1, expand=True)[0].values[vid_mask]

    return vids, rbp_mask, np.array(labels)[rbp_mask].tolist()


def get_plof_id_func(vep_tsv):
    # returns a function that takes an array of variant ids as input and returns a boolean array which indicates whether a variant is
    # annotated as a protein LOF variant for any gene

    plof_ids = pd.read_csv(vep_tsv, sep='\t', usecols=['Uploaded_variation'], index_col='Uploaded_variation').index.drop_duplicates()

    def is_plof(vids):
        return np.array([vid in plof_ids for vid in vids])

    return is_plof


# start: kernel parameters / functions
def get_cholesky(S):
    try:
        chol = cholesky(S)
        flag = 1
    except LinAlgError:
        try:
            np.fill_diagonal(S, S.diagonal() + 1e-6)  # sometimes this saves it...
            S /= np.max(S)
            chol = cholesky(S)
            flag = 0
        except LinAlgError:
            chol = np.eye(len(S))
            flag = -1
    return chol, flag


def get_simil_from_vep(V):
    # cosine similarity
    return cosine_similarity(V)


def get_simil_from_pos(pos):
    # 0.5 at 50bp distance
    gamma = -1*np.log(0.5)/(50**2)
    return rbf_kernel(pos[:, np.newaxis], gamma=gamma)


def get_weights_from_vep(V):
    # max norm
    return np.sqrt(np.max(np.abs(V), axis=1))


def get_regions():
    # load the results, keep those below a certain p-value
    results = pd.read_csv(snakemake.input.results_tsv, sep='\t', comment='#')

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


# end: kernel parameters / functions

# genotype path, vep-path:
assert len(snakemake.params.ids) == len(snakemake.input.bed), 'Error: length of chromosome IDs does not match length of genotype files'
geno_vep = zip(snakemake.params.ids, snakemake.input.bed, snakemake.input.mac_report, snakemake.input.ensembl_vep_tsv)

# get the top hits
regions_all = get_regions()
if regions_all is None:
    logging.info('No genes pass significance threshold, exiting.')
    sys.exit(0)

# where we store the results
stats = []
i_gene = 0

# dict with paths to input files
vep_h5 = {chrom: {'plus': p, 'minus': m} for chrom, (p, m) in zip(snakemake.params.ids, zip(snakemake.input.h5_rbp_plus, snakemake.input.h5_rbp_minus))}
vep_bed = {chrom: {'plus': p, 'minus': m} for chrom, (p, m) in zip(snakemake.params.ids, zip(snakemake.input.bed_rbp_plus, snakemake.input.bed_rbp_minus))}


# enter the chromosome loop:
timer = Timer()
for i, (chromosome, bed, mac_report, vep_tsv) in enumerate(geno_vep):

    if chromosome.replace('chr','') not in regions_all.chrom.unique():
        continue

    # get variants that pass MAF threshold:
    mac_report = maf_filter(mac_report)
    filter_vids = mac_report.index.values

    # function to identify protein LOF variants
    is_plof = get_plof_id_func(vep_tsv)

    # set up local collapsing
    # collapser = LocalCollapsing(distance_threshold=51.)

    for strand in ['plus', 'minus']:

        # set up the regions to loop over for the chromosome
        chromosome_id = chromosome.replace('chr', '')

        regions = regions_all[(regions_all.chrom == chromosome_id) & (regions_all.strand == strand)]

        if len(regions) == 0:
            continue

        # get variants that pass variant effect prediction threshold:
        vep_vids, vep_mask, labels = vep_filter(vep_h5[chromosome][strand], vep_bed[chromosome][strand])

        # combine
        filter_vids_chromosome = intersect_ids(vep_vids, filter_vids)

        # initialize the vep loader
        veploader = Hdf5Loader(vep_bed[chromosome][strand], vep_h5[chromosome][strand], 'diffscore', from_janggu=True)
        veploader.update_variants(filter_vids_chromosome)
        veploader.set_mask(vep_mask)

        # set up the variant loader (rbp variants) for the chromosome + strand
        plinkloader = VariantLoaderSnpReader(Bed(bed, count_A1=True, num_threads=4))
        plinkloader.update_variants(veploader.get_vids())
        plinkloader.update_individuals(covariatesloader.get_iids())


        # set up the genotype + vep loading function
        def get_rbp(interval):

            temp_genotypes, temp_vids, pos = plinkloader.genotypes_by_region(interval, return_pos=True)

            if temp_genotypes is None:
                raise GotNone

            ncarrier = np.nansum(temp_genotypes, axis=0)
            ncarrier = np.minimum(ncarrier, temp_genotypes.shape[0] - ncarrier).astype(int)

            temp_genotypes -= np.nanmean(temp_genotypes, axis=0)
            G1 = np.ma.masked_invalid(temp_genotypes).filled(0.)

            # deepripe variant effect predictions (single RBP)
            V1 = veploader.anno_by_id(temp_vids)

            weights = get_weights_from_vep(V1)

            S = get_simil_from_pos(pos)
            S *= get_simil_from_vep(V1) # Schur product of two positive definite matrices is positive definite

            cummac = mac_report.loc[temp_vids].Minor

            V1 = pd.DataFrame(V1, columns=labels)[sorted(labels)]

            return G1, temp_vids, weights, S, ncarrier, cummac, pos, V1


        # set up the test-function for a single gene
        def test_gene(interval, seed):

            interval = interval.to_dict()

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
                pval_dict[name]['pv_score'] = np.array([pv_score(GV[:, i, np.newaxis]) for i in range(GV.shape[1])])
                
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
            G, vids, weights, S, ncarrier, cummac, pos, V = get_rbp(interval)
            keep = ~is_plof(vids)

            # these are common to all kernels
            pval_dict['vid'] = vids
            pval_dict['weights'] = weights ** 2 # get_weights_from_vep returns the square root!
            pval_dict['MAC'] = cummac
            pval_dict['nCarrier'] = ncarrier
            pval_dict['not_LOF'] = keep
            for col in V.columns:
                pval_dict[col] = V[col].values.astype(np.float32)

            # # cholesky
            # if G.shape[1] > 1:
            #     L, flag1 = get_cholesky(S)
            # else:
            #     L, flag1 = np.eye(G.shape[1]), -1
            #
            # # do a score test (cholesky, and weighted cholesky)
            # GL = G.dot(L)
            # call_score(GL, 'lincholesky')
            # GWL = G.dot(np.diag(weights, k=0)).dot(L)
            # call_score(GWL, 'linwcholesky')

            # single-variant p-values:
            call_score(G, 'linw', vids=vids) # get's the single-variant p-values and effect-sizes -> no weighting!
            call_lrt(G.dot(np.diag(weights, k=0)), 'linw') # performs weighting and joint estimation of effect-sizes 

            # sanity checks
            assert len(vids) == interval['n_snp'], 'Error: number of variants does not match! expected: {}  got: {}'.format(interval['n_snp'], len(vids))
            assert cummac.sum() == interval['cumMAC'], 'Error: cumMAC does not match! expeced: {}, got: {}'.format(interval['cumMAC'], cummac.sum())

            if np.any(keep):

                # if keep.sum() == 1:
                #     # only single SNP is not LOF
                #     GL = G[:, keep] # actually just the linear kernel
                #     GWL = G[:, keep].dot(np.diag(weights[keep], k=0)) # actually just the linear weighted kernel
                # else:
                #     L, flag2 = get_cholesky(S[np.ix_(keep, keep)])
                #     GL = G[:, keep].dot(L)
                #     GWL = G[:, keep].dot(np.diag(weights[keep], k=0)).dot(L)

                # call_score(G1, 'lincholesky_notLOF')
                # call_score(GWL, 'linwcholesky_notLOF')
                #
                # call_lrt(GL, 'lincholesky_notLOF')
                # call_lrt(GWL, 'linwcholesky_notLOF')

                G1 = G[:,keep].dot(np.diag(weights[keep], k=0))
                call_score(G1, 'linw_notLOF', vids=vids[keep])
                call_lrt(G1, 'linw_notLOF')

            return pval_dict, called


        logging.info('loaders for chromosome {}, strand "{}" initialized in {:.1f} seconds.'.format(chromosome, strand, timer.check()))
        # run tests for all genes on the chromosome / strand
        for _, region in regions.iterrows():

            try:
                gene_stats, called = test_gene(region, i_gene)
            except GotNone:
                continue

            # build the single-variant datafame
            single_var_columns = ['gene', 'vid', 'weights', 'MAC', 'nCarrier', 'not_LOF' ] + sorted(snakemake.params.rbp_of_interest)
            sv_df = pd.DataFrame.from_dict({k: gene_stats[k] for k in single_var_columns})

            sv_df['pv_score'] = gene_stats['linw']['pv_score'] # single-variant p-values estimated independently
            sv_df['coef_random'] = gene_stats['linw']['coef_random'] # single-variant coefficients estimated jointly after weighting
            sv_df['beta'] = gene_stats['linw']['beta'] # single-variant coeffcients estimated independently *without* weighting
            sv_df['betaSd'] = gene_stats['linw']['betaSd'] # standard errors for the single-variant coefficients estimated independently *without* weighting
            sv_df['pheno'] = snakemake.params.phenotype

            out_dir = os.path.join(snakemake.params.out_dir_stats, region['name'])
            os.makedirs(out_dir, exist_ok=True)

            sv_df.to_csv(out_dir + '/variants.tsv.gz', sep='\t', index=False)

            for k in called:

                if k == 'variant_pvals':
                    continue

                results_dict = gene_stats[k]

                df_cols = ['pv_score', 'coef_random', 'beta', 'betaSd', 'vid']  # parts of the dict that have lenght > 1

                df = pd.DataFrame.from_dict(data={k: results_dict[k] for k in df_cols if k in results_dict})
                df['gene'] = gene_stats['gene']
                df['pheno'] = snakemake.params.phenotype
                df.to_csv(out_dir + '/{}.tsv.gz'.format(k), sep='\t', index=False)

                #  other cols ['nLL', 'sigma2', 'lrtstat', 'h2', 'log_delta']
                other_cols = {k: v for k, v in results_dict.items() if k not in df_cols}
                other_cols['gene'] = gene_stats['gene']
                other_cols['pheno'] = snakemake.params.phenotype

                pickle.dump(other_cols, open(out_dir + '/{}_stats.pkl'.format(k), 'wb'))

            i_gene += 1
            logging.info('tested {} genes...'.format(i_gene))

        timer.reset()
