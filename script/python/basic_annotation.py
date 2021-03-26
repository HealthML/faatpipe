#!/usr/bin/env python
# coding: utf-8

import pyranges
import pysam
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from pybedtools import BedTool, Interval
import os


def vcf_to_pyranges(vcf, tmpfile):
    '''
    create a bed-file containing the variant locations and ids
    '''
    
    bedtool = BedTool((Interval(record.chrom, record.pos - 1, record.pos - 1 + len(record.ref), name=record.id) for record in vcf.fetch() if not len(record.alts) > 1))
    bedtool.moveto(tmpfile)
    
    try:
        pr = pyranges.read_bed(tmpfile)
        del bedtool
        os.remove(tmpfile)
    except Exception as e:
        os.remove(tmpfile)
        raise e
    
    return pr


annot = pyranges.read_gtf(snakemake.input.anno_gtf)
annot = annot[annot.gene_biotype == 'protein_coding']
introns = annot.features.introns()

features = {k: annot[annot.Feature == k].drop() for k in ['CDS', 'five_prime_utr', 'three_prime_utr']}
features['intron'] = introns.drop()

# high confidence regions
highc = pyranges.read_bed(snakemake.input.hc_bed)



# exome sequencing target regions
exometarg = pyranges.read_bed(snakemake.input.es_bed)


# load the variants into a pyranges object
prpm = vcf_to_pyranges(pysam.VariantFile(snakemake.input.vcf), tmpfile = snakemake.output.tsv + '_tmp.bed')
# count overlaps to different feature types
prpm = pyranges.count_overlaps(features, prpm)

# annotation by majority vote
main_anno = np.argmax(prpm.as_df()[['CDS', 'five_prime_utr', 'three_prime_utr', 'intron']].values, axis=1)
d = {i: k for i, k in enumerate(['CDS', 'five_prime_utr', 'three_prime_utr', 'intron'])}
main_anno = pd.Series(main_anno).map(d)
prpm.major_anno = main_anno

# annotate as intron only those regions that have no other annotations
main_anno_2 = np.argmax( prpm.as_df()[['CDS', 'five_prime_utr', 'three_prime_utr']].values, axis=1)
i_anno = (prpm.as_df()[['CDS', 'five_prime_utr', 'three_prime_utr', 'intron']] > 0.).astype(int)
prpm.anno = pd.Series(np.where(np.sum(i_anno.iloc[:, 0:3], axis=1) >= i_anno['intron'], main_anno_2, 3)).map(d)

# overlap with exome target sequences
prpm = prpm.count_overlaps(exometarg, overlap_col='exomeseq_target')
prpm.exomeseq_target = (prpm.exomeseq_target > 0.).astype(int)

assert np.sum(prpm.exomeseq_target) > 0, 'Error: no variants found in exome-sequencing target regions. This could be a problem with incompatible chromosome names.'

# overlap with high-confidence genotyping regions
prpm = prpm.count_overlaps(highc, overlap_col='hiconf_reg')
prpm.hiconf_reg = (prpm.hiconf_reg > 0.).astype(int)

assert np.sum(prpm.hiconf_reg) > 0, 'Error: no variants found in high-confidence regions. This could be a problem with incompatible chromosome names.'


prpm.as_df()[['Name', 'major_anno', 'anno', 'CDS', 'five_prime_utr', 'three_prime_utr', 'intron', 'exomeseq_target', 'hiconf_reg']].to_csv(snakemake.output.tsv, sep='\t', index=False, header=True)

