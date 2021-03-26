
import os
# os.environ["OMP_NUM_THREADS"] = "16"

import logging
logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)

import pandas as pd
import numpy as np
import h5py

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

    rbp_mask = np.array([ i in rbp_of_interest for i in labels ])

    if np.ndim(diffscores) == 1:
        diffscores = diffscores[:, np.newaxis]

    vid_mask = np.max(np.abs(diffscores), axis=1) >= snakemake.params.min_impact

    vid_df = pd.read_csv(bed_rbp, sep='\t', header=None, usecols=[3], names=['vid'], dtype={'vid': str})
    vids = vid_df.vid.str.split(pat='_[ACGT]+>[ACGT]+$', n=1, expand=True)[0].values[vid_mask]

    return vids, rbp_mask


def get_plof_id_func(vep_tsv):

    # returns a function that takes an array of variant ids as input and returns a boolean array which indicates whether a variant is
    # annotated as a protein LOF variant for any gene

    plof_ids = pd.read_csv(vep_tsv, sep='\t', usecols=['Uploaded_variation'], index_col='Uploaded_variation').index.drop_duplicates()

    def is_plof(vids):
        return np.array([vid in plof_ids for vid in vids])

    return is_plof


# genotype path, vep-path:
assert len(snakemake.params.ids) == len(snakemake.input.bed), 'Error: length of chromosome IDs does not match length of genotype files'

geno_vep = zip(snakemake.params.ids, snakemake.input.bed, snakemake.input.mac_report, snakemake.input.ensembl_vep_tsv)

stats = []
simulations = []

i_gene = 0

# dict with paths to input files
vep_h5 = { chrom: {'plus': p, 'minus': m } for chrom, (p, m) in zip(snakemake.params.ids, zip(snakemake.input.h5_rbp_plus, snakemake.input.h5_rbp_minus))}
vep_bed = { chrom: {'plus': p, 'minus': m } for chrom, (p, m) in zip(snakemake.params.ids, zip(snakemake.input.bed_rbp_plus, snakemake.input.bed_rbp_minus))}

# regions
regions_all = pd.read_csv(snakemake.input.regions_bed, sep='\t', header=None, usecols=[0 ,1 ,2 ,3, 5], dtype={0: str, 1: np.int32, 2: np.int32, 3: str, 5: str})
regions_all.columns = ['chrom', 'start', 'end', 'name', 'strand']
regions_all['strand'] = regions_all.strand.map({'+': 'plus', '-': 'minus'})

# enter the chromosome loop:
timer = Timer()
for i, (chromosome, bed, mac_report, vep_tsv) in enumerate(geno_vep):

    if snakemake.params.debug:
        # skip most chromosomes if we are debugging...
        if chromosome not in ['chr9','chr16','chr21']:
            continue

    # get variants that pass MAF threshold:
    mac_report = maf_filter(mac_report)
    filter_vids = mac_report.index.values

    # function to identify protein LOF variants
    is_plof = get_plof_id_func(vep_tsv)

    # set up local collapsing
    collapser = LocalCollapsing(distance_threshold=51.)

    for strand in ['plus', 'minus']:

        # set up the regions to loop over for the chromosome
        chromosome_id = chromosome.replace('chr','')

        regions = regions_all[(regions_all.chrom == chromosome_id) & (regions_all.strand == strand)]

        # get variants that pass variant effect prediction threshold:
        vep_vids, vep_mask = vep_filter(vep_h5[chromosome][strand], vep_bed[chromosome][strand])
        
        # assert len(vep_vids) > 0

        # combine
        filter_vids_chromosome = intersect_ids(vep_vids, filter_vids)
        
        # assert len(filter_vids_chromosome) > 0

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
            weights = V1[:,0]

            cummac = mac_report.loc[temp_vids].Minor

            return G1, temp_vids, weights, ncarrier, cummac, pos


        # set up the test-function for a single gene
        def test_gene(interval, seed):
            
            interval = interval.to_dict()

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
            G, vids, weights, ncarrier, cummac, pos = get_rbp(interval)
            keep = ~is_plof(vids)

            # perform local collapsing with weights
            if G.shape[1] > 1:
                G1, clusters1 = collapser.collapse(G, pos, np.sign(weights) * np.sqrt(np.abs(weights)))
            else:
                G1 = G.dot(np.diag(np.sqrt(np.abs(weights)), k=0))
                clusters1 = [0]

            # do a score test (local collapsing)
            call_score(G1, 'linwcollapsed')

            clusters2 = None
            # if gene is nominally significant:
            if pval_dict['pv_score_linwcollapsed'] < 0.01:

                # do lrt tests
                call_lrt(G1, 'linwcollapsed')

                if np.any(keep):

                    if keep.sum() == 1:
                        # only single SNP is not LOF
                        G2 = G[:, keep].dot(np.diag(np.sqrt(np.abs(weights[keep])), k=0))
                        clusters2 = [0]
                    else:
                        G2, clusters2 = collapser.collapse(G[:,keep], pos[keep], np.sign(weights[keep]) * np.sqrt(np.abs(weights[keep])))

                    call_score(G2, 'linwcollapsed_notLOF')
                    call_lrt(G2, 'linwcollapsed_notLOF')


            pval_dict['nCarrier'] = ncarrier.sum()
            pval_dict['cumMAC'] = cummac.sum()
            pval_dict['n_snp'] = len(vids)
            pval_dict['n_cluster'] = len(set(clusters1))

            pval_dict['n_snp_notLOF'] = keep.sum()
            pval_dict['n_cluster_notLOF'] = len(set(clusters2)) if clusters2 is not None else -1

            # sanity check
            assert pval_dict['cumMAC'] >= pval_dict['nCarrier'], 'Error: something is broken.'

            return pval_dict, sim_dict

        logging.info('loaders for chromosome {}, strand "{}" initialized in {:.1f} seconds.'.format(chromosome, strand, timer.check()))
        # run tests for all genes on the chromosome / strand
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

    logging.info('all tests for chromosome {} performed in {:.2f} minutes.'.format(chromosome, timer.check( ) /60.))
    logging.info('genes tested so far: {}'.format(i_gene))


# when all chromosomes are done:

# generate results table:
stats = pd.DataFrame.from_dict(stats).set_index('gene')

# these are the kernels
kern = ['linwcollapsed', 'linwcollapsed_notLOF']

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
pvals['rbp'] = snakemake.params.rbp_of_interest[0]

pvals.to_csv(snakemake.output.results_tsv, sep='\t', index=True)

params = pd.DataFrame(params)
params.to_csv(snakemake.output.chi2param_tsv, sep='\t', index=True)


