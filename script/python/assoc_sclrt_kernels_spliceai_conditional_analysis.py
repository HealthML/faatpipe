import os
# os.environ["OMP_NUM_THREADS"] = "16"

import logging

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO)

import pandas as pd
import numpy as np

# seak imports
from seak.data_loaders import intersect_ids, EnsemblVEPLoader, VariantLoaderSnpReader, CovariatesLoaderCSV
from seak.scoretest import ScoretestNoK
from seak.lrt import LRTnoK
import gzip

from pysnptools.snpreader import Bed
import pickle
import sys

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

# conditional analysis:
# we need to fit gene-specific null models for the LRT!
# null_model_lrt = LRTnoK(X, Y)
null_model_lrt = None


# set up function to filter variants:
def maf_filter(mac_report):
    # load the MAC report, keep only observed variants with MAF below threshold
    mac_report = pd.read_csv(mac_report, sep='\t', usecols=['SNP', 'MAF', 'Minor', 'alt_greater_ref'])

    if snakemake.params.filter_highconfidence:
        vids = mac_report.SNP[(mac_report.MAF < snakemake.params.max_maf) & (mac_report.Minor > 0) & ~(mac_report.alt_greater_ref.astype(bool)) & (mac_report.hiconf_reg.astype(bool))]
    else:
        vids = mac_report.SNP[(mac_report.MAF < snakemake.params.max_maf) & (mac_report.Minor > 0) & ~(mac_report.alt_greater_ref.astype(bool))]

    # this has already been done in filter_variants.py
    # load the variant annotation, keep only variants in high-confidece regions
    # anno = pd.read_csv(anno_tsv, sep='\t', usecols=['Name', 'hiconf_reg'])
    # vids_highconf = anno.Name[anno.hiconf_reg.astype(bool).values]
    # vids = np.intersect1d(vids, vids_highconf)

    return mac_report.set_index('SNP').loc[vids]

def sid_filter(vids):
    
    if 'sid_include' in snakemake.config:
        print('limiting to variants present in {}'.format(snakemake.config['sid_include']))
        
        infilepath = snakemake.config['sid_include']
        
        if infilepath.endswith('gz'):
            with gzip.open(infilepath,'rt') as infile:
                sid = np.array([l.rstrip() for l in infile])
        else:
            with open(infilepath, 'r') as infile:
                sid = np.array([l.rstrip() for l in infile])
    else:
        return vids
                
    return intersect_ids(vids, sid)


def get_regions():
    # load the results, keep those below a certain p-value + present in conditional analysis list 
    results = pd.read_csv(snakemake.input.results_tsv, sep='\t')

    kern = snakemake.params.kernels
    if isinstance(kern, str):
        kern = [kern]

    pvcols_score = ['pv_score_' + k for k in kern]
    pvcols_lrt = ['pv_lrt_' + k for k in kern]
    statcols = ['lrtstat_' + k for k in kern]
    results = results[['gene', 'n_snp', 'cumMAC', 'nCarrier'] + statcols + pvcols_score + pvcols_lrt]
    
    results.rename(columns={'gene':'name'}, inplace=True)
    
    # set up the regions to loop over for the chromosome
    regions = pd.read_csv(snakemake.input.regions_bed, sep='\t', header=None, usecols=[0 ,1 ,2 ,3, 5], dtype={0 :str, 1: np.int32, 2 :np.int32, 3 :str, 5:str})
    regions.columns = ['chrom', 'start', 'end', 'name', 'strand']
    regions['strand'] = regions.strand.map({'+': 'plus', '-': 'minus'})
    
    # conditional analysis, keep only genes that are provided in conditional_list
    regions['gene_name'] = regions['name'].str.split('_', expand=True)[1]
    
    _conditional_list = pd.read_csv(snakemake.input.conditional_list, sep='\t', header=0, usecols=['gene_name','pheno'])
    _conditional_list = _conditional_list[(_conditional_list.pheno == snakemake.wildcards['pheno']) | (_conditional_list.pheno == snakemake.params.phenotype) ]
    _conditional_list = _conditional_list.drop_duplicates() 
    
    if len(_conditional_list) == 0:
        logging.info('No genes pass significance threshold for phenotype {}, exiting.'.format(snakemake.params.phenotype))
        with gzip.open(snakemake.output.results_tsv, 'wt') as outfile:
            outfile.write('# no significant hits for phenotype {}. \n'.format(snakemake.params.phenotype))
        sys.exit(0)
    else:
        regions = regions.merge(_conditional_list, how='right')
        
    regions = regions.merge(results, how='left', on=['name'], validate='one_to_one')
    
    sig_genes = [results.name[results[k] < snakemake.params.significance_cutoff].values for k in pvcols_score + pvcols_lrt]
    sig_genes = np.unique(np.concatenate(sig_genes))
    
    if len(sig_genes) == 0:
        return None
    
    regions = regions.loc[regions.name.isin(sig_genes)]
        
    return regions

# genotype path, vep-path:
assert len(snakemake.params.ids) == len (snakemake.input.bed), 'Error: length of chromosome IDs does not match length of genotype files'
geno_vep = zip(snakemake.params.ids, snakemake.input.bed, snakemake.input.vep_tsv, snakemake.input.ensembl_vep_tsv, snakemake.input.mac_report, snakemake.input.h5_lof, snakemake.input.iid_lof, snakemake.input.gid_lof)

# get the top hits
regions_all = get_regions()
if regions_all is None:
    logging.info('No genes pass significance threshold, exiting.')
    import gzip
    with gzip.open(snakemake.output.results_tsv, 'wt') as outfile:
        outfile.write('# no hits below specified threshold ({}) for phenotype {}. \n'.format(snakemake.params.significance_cutoff, snakemake.params.phenotype))
    sys.exit(0)

logging.info('About to evaluate variants in {} genes'.format(len(regions_all.name.unique())))

# storing all results here
results = []

i_gene = 0
i_chrom = 0

# conditional analysis:
# these genotypes will be used for the conditional analysis
geno_cond = VariantLoaderSnpReader(Bed(snakemake.input.conditional_geno, count_A1=True, num_threads=1))
geno_cond.update_individuals(covariatesloader.get_iids())

# this file contains the mapping of associations to SNPs to condition on
conditional_list = pd.read_csv(snakemake.input.conditional_list, sep='\t', header=0)
conditional_list = conditional_list[(conditional_list.pheno == snakemake.params['phenotype']) | (conditional_list.pheno == snakemake.wildcards['pheno']) ].drop_duplicates()

geno_cond.update_variants(intersect_ids(conditional_list.snp_rsid, geno_cond.get_vids()))
logging.info('considering {} variants as covariates for conditional tests.'.format(len(geno_cond.get_vids())))

# enter the chromosome loop:
timer = Timer()
for i, (chromosome, bed, vep_tsv, ensembl_vep_tsv, mac_report, h5_lof, iid_lof, gid_lof) in enumerate(geno_vep):

    if chromosome.replace('chr','') not in regions_all.chrom.unique():
        continue

    if snakemake.params.debug:
        # process only two chromosomes if debugging...
        if i_chrom > 2:
            break

    timer.reset()

    # set up the ensembl vep loader for the chromosome
    spliceaidf = pd.read_csv(vep_tsv,
                             sep='\t',
                             usecols=['name', 'chrom', 'end', 'gene', 'max_effect', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL'],
                             index_col='name')

    # get set of variants for the chromosome:
    mac_report = maf_filter(mac_report)
    filter_vids = mac_report.index.values
    filter_vids = sid_filter(filter_vids)

    # filter by MAF
    keep = intersect_ids(filter_vids, spliceaidf.index.values)
    spliceaidf = spliceaidf.loc[keep]
    spliceaidf.reset_index(inplace=True)

    # filter by impact:
    spliceaidf = spliceaidf[spliceaidf.max_effect >= snakemake.params.min_impact]

    # set up the regions to loop over for the chromosome
    regions = regions_all.copy()

    # discard all genes for which we don't have annotations
    gene_ids = regions.name.str.split('_', expand=True)  # table with two columns, ensembl-id and gene-name
    regions['gene'] = gene_ids[1]  # this is the gene name
    regions['ensembl_id'] = gene_ids[0]
    regions.set_index('gene', inplace=True)
    genes = intersect_ids(np.unique(regions.index.values), np.unique(spliceaidf.gene))  # intersection of gene names
    regions = regions.loc[genes].reset_index()  # subsetting
    regions = regions.sort_values(['chrom', 'start', 'end'])

    # check if the variants are protein LOF variants, load the protein LOF variants:
    ensemblvepdf = pd.read_csv(ensembl_vep_tsv, sep='\t', usecols=['Uploaded_variation', 'Gene'])

    # this column will contain the gene names:
    genes = intersect_ids(np.unique(ensemblvepdf.Gene.values), regions.ensembl_id)  # intersection of ensembl gene ids
    ensemblvepdf = ensemblvepdf.set_index('Gene').loc[genes].reset_index()
    ensemblvepdf['gene'] = gene_ids.set_index(0).loc[ensemblvepdf.Gene.values].values

    # set up the merge
    ensemblvepdf.drop(columns=['Gene'], inplace=True)  # get rid of the ensembl ids, will use gene names instead
    ensemblvepdf.rename(columns={'Uploaded_variation': 'name'}, inplace=True)
    ensemblvepdf['is_plof'] = 1.

    ensemblvepdf = ensemblvepdf[~ensemblvepdf.duplicated()]  # if multiple ensembl gene ids map to the same gene names, this prevents a crash.

    # we add a column to the dataframe indicating whether the variant is already annotated as protein loss of function by the ensembl variant effect predictor
    spliceaidf = pd.merge(spliceaidf, ensemblvepdf, on=['name', 'gene'], how='left', validate='one_to_one')
    spliceaidf['is_plof'] = spliceaidf['is_plof'].fillna(0.).astype(bool)

    # initialize the loader
    # Note: we use "end" here because the start + 1 = end, and we need 1-based coordiantes (this would break if we had indels)
    eveploader = EnsemblVEPLoader(spliceaidf['name'], spliceaidf['chrom'].astype('str') + ':' + spliceaidf['end'].astype('str'), spliceaidf['gene'], data=spliceaidf[['max_effect', 'is_plof', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL']].values)

    # set up the variant loader (splice variants) for the chromosome
    plinkloader = VariantLoaderSnpReader(Bed(bed, count_A1=True, num_threads=4))
    plinkloader.update_variants(eveploader.get_vids())
    plinkloader.update_individuals(covariatesloader.get_iids())

    # set up the protein LOF burden loader
    bloader_lof = BurdenLoaderHDF5(h5_lof, iid_lof, gid_lof)
    bloader_lof.update_individuals(covariatesloader.get_iids())


    # set up the splice genotype + vep loading function
    def get_splice(interval):

        try:
            V1 = eveploader.anno_by_interval(interval, gene=interval['name'].split('_')[1])
        except KeyError:
            raise GotNone

        if V1.index.empty:
            raise GotNone

        vids = V1.index.get_level_values('vid')
        V1 = V1.droplevel(['gene'])

        temp_genotypes, temp_vids = plinkloader.genotypes_by_id(vids, return_pos=False)

        ncarrier = np.nansum(temp_genotypes > 0, axis=0)

        temp_genotypes -= np.nanmean(temp_genotypes, axis=0)
        G1 = np.ma.masked_invalid(temp_genotypes).filled(0.)

        cummac = mac_report.loc[vids].Minor

        # spliceAI max score
        weights = V1[0].values.astype(np.float64)

        is_plof = V1[1].values.astype(bool)

        splice_preds_all = V1.iloc[:,2:]
        splice_preds_all.columns = ['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL']

        # "standardized" positions -> codon start positions
        # pos = V1[0].values.astype(np.int32)

        return G1, vids, weights, ncarrier, cummac, is_plof, splice_preds_all


    # set up the protein-LOF loading function

    def get_plof(interval):

        try:
            G2 = bloader_lof.genotypes_by_id(interval['name']).astype(np.float64)
        except KeyError:
            G2 = None

        return G2

    
    # conditional analysis
    # set up the function to load the variants to condition on
    
    def get_conditional(interval):
        
        cond_snps_vid = conditional_list.loc[ conditional_list.gene_name == interval['gene_name'], 'snp_rsid' ].unique().tolist()
        
        temp_genotypes, temp_vids = geno_cond.genotypes_by_id(cond_snps_vid, return_pos=False)
        temp_genotypes -= np.nanmean(temp_genotypes, axis=0)
        
        cond_snps = np.ma.masked_invalid(temp_genotypes).filled(0.)

        return cond_snps, temp_vids
    

    # set up the test-function for a single gene
    def test_gene(interval, seed):

        pval_dict = {}
        pval_dict['gene'] = interval['name']

        out_dir = os.path.join(snakemake.params.out_dir_stats, interval['name'])
        os.makedirs(out_dir, exist_ok=True)

        # conditional analysis:
        # get the snps to condition on, and include them in the null model for the LRT
        cond_snps, cond_snps_vid = get_conditional(interval)
        null_model_lrt = LRTnoK(np.concatenate([X, cond_snps], axis=1), Y)
        
        # conditional analysis:
        # the score-test takes a second argument (G2) that allows conditioning on a second set of variants...
        def pv_score(GV, G2=cond_snps):
            # wraps score-test
            pv = null_model_score.pv_alt_model(GV, G2)
            if pv < 0.:
                pv = null_model_score.pv_alt_model(GV, G2, method='saddle')
            return pv

        def call_test(GV, name):
            pval_dict['pv_score_' + name] = pv_score(GV)
            
            altmodel = null_model_lrt.altmodel(GV)

            res = null_model_lrt.pv_sim_chi2(250000, simzero=False, seed=seed)
            pval_dict['pv_lrt_' + name] = res['pv']
            pval_dict['lrtstat_' + name ] = altmodel['stat']
            if 'h2' in altmodel:
                pval_dict['h2_' + name ] = altmodel['h2']
            if res['pv'] != 1.:
                for stat in ['scale', 'dof', 'mixture', 'imax']:
                    pval_dict[stat + '_' + name] = res[stat]
                if len(res['res'] > 0):
                    pd.DataFrame({interval['name']: res['res']}).to_pickle(out_dir + '/{}.pkl.gz'.format(name))


        # load splice variants
        G1, vids, weights, ncarrier, cummac, is_plof, splice_preds_all = get_splice(interval)
        # keep indicates which variants are NOT "protein LOF" variants, i.e. variants already identified by the ensembl VEP
        keep = ~is_plof

        # sanity checks
        assert len(vids) == interval['n_snp'], 'Error: number of variants does not match! expected: {}  got: {}'.format(interval['n_snp'], len(vids))
        assert cummac.sum() == interval['cumMAC'], 'Error: cumMAC does not match! expeced: {}, got: {}'.format(interval['cumMAC'], cummac.sum())


        # do a score burden test (max weighted), this is different than the baseline!
        G1_burden = np.max(np.where(G1 > 0.5, np.sqrt(weights), 0.), axis=1, keepdims=True)
        
        try:
            call_test(G1_burden, 'linwb')
        except AssertionError as err:
            out_dump_np = '{}_{}_{{}}_covar.npy'.format(interval['name'],snakemake.params.phenotype)
            np.save(out_dump_np.format('G1burden'), G1_burden)
            np.save(out_dump_np.format('G1'), G1)
            np.save(out_dump_np.format('Covar'), np.concatenate([X, cond_snps], axis=1))
            np.save(out_dump_np.format('Y'), Y)
            logging.error('AssertionError encountered when testing gene {}, conditioning on {}. dumping Y, covariates and G1 to {}'.format(interval['name'], ','.join(cond_snps_vid),out_dump_np.format('*')))
            raise err

        # linear weighted kernel
        G1 = G1.dot(np.diag(np.sqrt(weights), k=0))

        # do a score test (linear weighted)
        call_test(G1, 'linw')

        # load plof burden
        G2 = get_plof(interval)

        if G2 is not None:

            if np.any(keep):

                # merged (single variable)
                G1_burden_mrg = np.maximum(G2, G1_burden)
                call_test(G1_burden_mrg, 'linwb_mrgLOF')

                # concatenated ( >= 2 variables)
                # we separate out the ones that are already part of the protein LOF variants!

                G1 = np.concatenate([G1[:, keep], G2], axis=1)
                call_test(G1, 'linw_cLOF')
            else:
                logging.info('All Splice-AI variants for gene {} where already identified by the Ensembl variant effect predictor'.format(interval['name']))

        # conditional analysis: keep names of SNPs that we condition on 
        pval_dict['cond_snps'] = ','.join(cond_snps_vid)
                
        return pval_dict


    logging.info('loaders for chromosome {} initialized in {:.1f} seconds.'.format(chromosome, timer.check()))
    # run tests for all genes on the chromosome
    for _, region in regions.iterrows():

        if snakemake.params.debug:
            if i_gene > 2:
                continue
        
        try:
            results.append(test_gene(region, i_gene)) # i_gene is used as the random seed to make things reproducible
        except GotNone:
            continue

        i_gene += 1
        logging.info('tested {} genes...'.format(i_gene))

    i_chrom += 1
    timer.reset()

# export the results in a single table
results = pd.DataFrame(results)
results['pheno'] = snakemake.params.phenotype
results.to_csv(snakemake.output.results_tsv, sep='\t', index=False)
