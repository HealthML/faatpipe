#!/usr/bin/env python

import allel
import numpy as np
import pandas as pd
from argparse import ArgumentParser


'''
deprecated.
'''

p = ArgumentParser()
p.add_argument('-vcf', help='input VCF file', required=True)
p.add_argument('-out', help='output TSV file', required=True)
p.add_argument('-chr', help='chromosome', required=True)
args=p.parse_args()

path_to_vcf_splice_scores = args.vcf
csv_path = args.out

chromosome = args.chr

if chromosome.startswith('chr'):
    chromosome = chromosome.replace('chr','')

spliceai_vcf_df = allel.vcf_to_dataframe(path_to_vcf_splice_scores,
                                         fields=['variants/CHROM',
                                                 'variants/POS',
                                                 'variants/ID',
                                                 'variants/REF',
                                                 'variants/ALT',
                                                 'variants/SYMBOL',
                                                 'variants/STRAND',
                                                 'variants/TYPE',
                                                 'variants/DIST',
                                                 'variants/DS_AG',
                                                 'variants/DS_AL',
                                                 'variants/DS_DG',
                                                 'variants/DS_DL',
                                                 'variants/DP_AG',
                                                 'variants/DP_AL',
                                                 'variants/DP_DG',
                                                 'variants/DP_DL',
                                                 'variants/is_snp'],
                                         alt_number=1,
                                         region=chromosome
                                         )

# Give every variant an ID that is the same as in the binary PLINK genotype files (bim files)
spliceai_vcf_df['ID'] = spliceai_vcf_df['CHROM'].astype(str) + ':' + spliceai_vcf_df['POS'].astype(str) + ':' \
                        + spliceai_vcf_df['REF'].astype(str) + ':' + spliceai_vcf_df['ALT'].astype(str)

spliceai_vcf_df['start'] = spliceai_vcf_df['POS'] - 1
spliceai_vcf_df['end'] = spliceai_vcf_df['POS']

# meta-info file with: 'chrom', 'start', 'end', 'name' (as in UCSC BED file) 0-based half-open; and additional columns
# 'gene' and 'strand'
# (2)

csv_path = csv_path + '.gz' if not csv_path.endswith('.gz') else csv_path

spliceai_vcf_df[['CHROM', 'start', 'end', 'ID', 'SYMBOL', 'STRAND', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL', 'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL']].rename(columns = {'CHROM':'#CHROM'}).to_csv(csv_path, header=True, index=False, sep='\t')



