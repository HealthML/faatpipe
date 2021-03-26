#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import logging

# filtering and overlapping the spliceAI variants with the UKB exome sequencing variants:

logging.basicConfig(filename=snakemake.log[0])

wgscores = pd.read_csv(snakemake.input.tsv, sep='\t', header=0, dtype={'chrom':str})


# reorder
wgscores = wgscores[['chrom', 'start', 'end', 'name', 'gene','strand', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']]

# get max effect
wgscores['max_effect']=np.max(wgscores[['DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']], axis=1)


# remove duplicates, keep only those with the maximum effect
wgscores = wgscores.loc[ wgscores.max_effect == wgscores.groupby(['name','gene'])['max_effect'].transform(np.max)]
wgscores = wgscores[~wgscores[['name','gene']].duplicated(keep=False)]

# load the variants id from uk-biobank (note: this of course assumes the variant names are chr:start:ref:alt)
variants = pd.read_csv(snakemake.input.bim, sep='\t', header=None, usecols=[1])[1]

logging.info('ready to intersect with {} variants from {}.'.format(len(variants), snakemake.input.bim))

wgscores.set_index('name', inplace=True)
# perform the intersection
variants_f = np.intersect1d(wgscores.index.values, variants)

if len(variants_f) == 0:
    logging.warn('no variants overlapped spliceAI variants!')

wgscores = wgscores.loc[variants_f]
wgscores.reset_index(inplace=True)

wgscores.to_csv(snakemake.output.tsv, sep='\t', index=False)





