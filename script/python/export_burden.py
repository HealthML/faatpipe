

from seak.data_loaders import EnsemblVEPLoader, VariantLoaderSnpReader, BEDRegionLoader, intersect_ids
from pysnptools.snpreader import Bed
import pandas as pd
import numpy as np
import sys
import h5py
import logging
import os

class GotNone(Exception):
    pass

'''
Exports "burdens" (indicator variables 0/1) which can be used for downstream analysis
'''

logging.basicConfig(filename=snakemake.log[0])

def initialize_h5(path, n_genes, n_individuals):
    '''
    initialize the results-file
    '''

    resultsfile = h5py.File(path, 'w-')
    resultsfile.create_dataset('G', maxshape=(n_genes, n_individuals), data=np.zeros((n_genes, n_individuals), dtype='i1'), dtype='i1', compression="gzip")
    resultsfile.close()


def maf_filter():
    '''
    determine which variants to keep based on MAF
    '''

    # load the MAC report, keep only variants that are observed and below missingness threshold
    mac_report = pd.read_csv(snakemake.input.mac_report, sep='\t', usecols=['SNP', 'Minor', 'MAF', 'MISSING', 'alt_greater_ref'])
    vids = mac_report.SNP[(mac_report.Minor > 0) & (mac_report.MISSING < snakemake.params.max_missing) & (mac_report.MAF < snakemake.params.max_maf) & ~(mac_report.alt_greater_ref.astype(bool))].values

    # load the variant annotation, keep only variants in high-confidece regions
    if snakemake.params.filter_highconfidence:
        
        anno = pd.read_csv(snakemake.input.anno_tsv, sep='\t', usecols=['Name', 'hiconf_reg'])
        vids_highconf = anno.Name[anno.hiconf_reg.astype(bool).values]
        vids = np.intersect1d(vids, vids_highconf)

    return vids

def iid_filter():
    '''
    determine which individuals to include
    '''

    # iid = pd.read_csv(snakemake.input.complete_cases, sep='\t', usecols=['iid'], header=0).iid.values()
    with open(snakemake.input.complete_cases, 'r') as infile:
        iid = [ l.rstrip() for l in infile ]
    return iid

def write_ids(ids, handle):
    for id in ids:
        handle.write('{}\n'.format(id))


def get_veploader_and_regions(filter_vids):

    '''
    get the ensembl-vep loader and gene regions
    '''

    if snakemake.params.effect == 'LOF':

        ensemblvepdf = pd.read_csv(snakemake.input.vep_tsv, sep='\t', usecols=['Uploaded_variation', 'Location', 'Gene'], index_col='Uploaded_variation')

        keep = intersect_ids(filter_vids, ensemblvepdf.index.values)
        ensemblvepdf = ensemblvepdf.loc[keep]
        ensemblvepdf.reset_index(inplace=True)

        eveploader = EnsemblVEPLoader(ensemblvepdf['Uploaded_variation'], ensemblvepdf['Location'], ensemblvepdf['Gene'])

    elif snakemake.params.effect == 'missense':

        ensemblvepdf = pd.read_csv(snakemake.input.vep_tsv, sep='\t', usecols=['Uploaded_variation', 'Location', 'Gene', 'pos_standardized', 'impact'], index_col='Uploaded_variation')

        keep = intersect_ids(filter_vids, ensemblvepdf.index.values)
        ensemblvepdf = ensemblvepdf.loc[keep]
        ensemblvepdf.reset_index(inplace=True)

        # filter by impact
        # since we can't weigh the variants. We use a different way to filter them than in the other tests:
        
        # NOT:
        # ensemblvepdf = ensemblvepdf[ensemblvepdf.groupby(['Gene', 'pos_standardized'])['impact'].transform(np.max) >= snakemake.params.min_impact]
        
        ensemblvepdf = ensemblvepdf[ ensemblvepdf.impact >= snakemake.params.min_impact ]
        eveploader = EnsemblVEPLoader(ensemblvepdf['Uploaded_variation'], ensemblvepdf['Location'], ensemblvepdf['Gene'], data=ensemblvepdf[['impact']].values)

    else:
        raise NotImplementedError('effect has to be either missense or LOF!')

    regions = pd.read_csv(snakemake.input.regions_bed, sep='\t', header=None, usecols=[0,1,2,3], dtype={0:str, 1: np.int32, 2:np.int32, 3:str})
    regions.columns = ['chrom', 'start', 'end', 'name']

    # discard all genes that are not on the chromosomes we are looking at:
    # chromosomes = np.unique(eveploader.pos_df.chrom)
    # regions = regions[regions.chrom.str.isin(chromosomes)]

    # discard all genes for which we don't have annotations
    regions['gene'] = regions.name.str.split('_', expand=True)[0]
    regions.set_index('gene', inplace=True)

    genes = intersect_ids(np.unique(regions.index.values), np.unique(eveploader.pos_df.gene))
    regions = regions.loc[genes].reset_index()

    regions = regions.sort_values(['chrom','start','end'])[['chrom','start','end','name','gene']]

    return eveploader, regions


def main():

    # get variants that pass MAF and genotyping filters
    filter_vids = maf_filter()

    # get the variant effect predictions and gene regions:
    eveploader, regions = get_veploader_and_regions(filter_vids)

    # get the genotype loader:
    plinkloader = VariantLoaderSnpReader(Bed(snakemake.input.genotypes_bed, count_A1=True))

    # intersect variants
    common_vids = intersect_ids(plinkloader.get_vids(), eveploader.get_vids())

    plinkloader.update_variants(common_vids)
    eveploader.update_variants(common_vids)

    # drop irrelevant indidivuals
    iids = iid_filter()
    plinkloader.update_individuals(iids)

    # batch size to write genotypes
    batch_size = 100

    def load_and_proc_geno(interval):

        try:
            V1 = eveploader.anno_by_interval(interval, gene=interval['name'].split('_')[0])
        except KeyError:
            raise GotNone

        if V1.index.empty:
            raise GotNone

        vids = V1.index.get_level_values('vid')

        temp_genotypes, temp_vids = plinkloader.genotypes_by_id(vids, return_pos=False)
        
        temp_genotypes = np.ma.masked_invalid(temp_genotypes).filled(0.) # since we have already kicked out all the "weird" variants, we can skip the pre-processing below

        #G1, vids = plinkloader.preprocess_genotypes(temp_genotypes,
        #                                              temp_vids,
        #                                              recode_maf=False,
        #                                              invert_encoding=False,
        #                                              impute_mean=True,
        #                                              center=True,
        #                                              max_maf=snakemake.params.max_maf)  # this will kick out any where major/minor are "flipped"

        #if G1 is None:
        #    raise GotNone

        G1_burden = (np.sum(temp_genotypes > 0.5, axis=1, keepdims=True) > 0.).astype('i1')

        return G1_burden, vids

    genos = []
    gene_ids = []

    # initialize output file
    initialize_h5(snakemake.output.h5, len(regions), len(iids))

    # the actual loop that glues it all together
    regions = regions.iterrows()

    h5 = h5py.File(snakemake.output.h5, 'r+')
    out = h5['G']
    
    idfile = open(snakemake.output.gene_txt, 'x')

    i = 0 # keeps track of how many entries have been exported to the hdf5 file
    ibatch = 0
    
    try:
        while True:

            b = []
            ids = []

            while ibatch < batch_size:

                _, region = next(regions)

                try:
                    G, vids = load_and_proc_geno(region)
                except GotNone:
                    continue

                b.append(G)
                ids.append(region['name'])

                ibatch += 1

            out[i:(i+ibatch)] = np.concatenate(b, axis = 1).T
            write_ids(ids, idfile)

            i += len(ids)
            ibatch = 0

    except StopIteration:

        if len(b) > 0:

            out[i:(i+ibatch)] = np.concatenate(b, axis = 1).T
            write_ids(ids, idfile)
            i += len(ids)

        idfile.close()
        out.resize(i, axis=0)
        h5.close()
    
    relpath_out = os.path.relpath(snakemake.input.complete_cases, os.path.dirname(snakemake.output.iid_txt) )
    
    os.symlink(relpath_out, snakemake.output.iid_txt)


if __name__ == '__main__':
    main()













