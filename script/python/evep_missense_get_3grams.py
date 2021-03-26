#!/usr/bin/env python
# coding: utf-8

# to be run in genomics environment

# extracts trimers of AA surrounding missense variants and computes the cosine similarity of the alt vs ref trimers

# input : ensembl variant effect predictor output

import pysam

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from argparse import ArgumentParser
from collections import defaultdict
import pandas as pd
import numpy as np
import sys
import gzip
from sklearn.metrics.pairwise import cosine_similarity

# https://github.com/remomomo/Deep-Proteomics.git
embed = pd.read_csv('./data/Deep-Proteomics/protVec_100d_3grams.csv', sep='\t', engine='python', header=None, index_col=0)

p = ArgumentParser()
p.add_argument('--infile', required=True)
p.add_argument('--out_prefix', required=True)
p.add_argument('--cds_fasta', required=True)
args = p.parse_args()

outfile = args.out_prefix + '.tsv' # save the results here
out_err = args.out_prefix + '_err.tsv' # save errors here
out_succ = args.out_prefix + '_succ.tsv' # save successes here

# cds fasta file 
cds = pysam.FastaFile(args.cds_fasta)

class Mismatch(Exception):
    pass

def get_coord(i, pseq):
    
    if i == 0:
        start = 0
        end = 3
        p = 0
    elif i == (len(pseq) - 2):
        if pseq.endswith('*'):
            start = i - 2
            end = i + 1
            p = 2
        else:
            start = i - 1
            end = i + 1
            p = 1
    elif i == (len(pseq) - 1):
        if pseq.endswith('*'):
            raise Mismatch
        start = i - 2
        end = i + 1
        p = 2
    else:
        start = i-1
        end = i+2
        p = 1

    return start, end, p


def zero():
    return 0

complement = {'A':'T','C':'G','G':'C','T':'A'}

def get_cosine(ref, alt):
    try:
        return cosine_similarity(embed.loc[[ref]],embed.loc[[alt]])[0][0]
    except KeyError:
        return 'nan'

def parse_scores(x):
    
    sift_str = ''
    sift_score = 'nan'
    polyphen_str = ''
    polyphen_score = 'nan'
    
    fields = x.split(';')
    
    for field in fields:
        if field.startswith('SIF'):
            sift_str, field_score = field[5:].split('(')
            sift_score = field_score[:-1]
        elif field.startswith('Pol'):
            polyphen_str, polyphen_score = field[9:].split('(')
            polyphen_score = polyphen_score[:-1]
            break
            
    return sift_str, sift_score, polyphen_str, polyphen_score
    

with open(outfile, 'x') as out:
    
    out.write('\t'.join(['Uploaded_variation','Location','Allele','Gene','Feature','Strand','Protein_position','codon_position','ref','alt','cosine_similarity','sift_str','sift_score','polyphen_str','polyphen_score\n']))
    
    with gzip.open(args.infile,'rt') as infile:
        i = 0
        
        g_succ = defaultdict(zero)
        g_err = defaultdict(zero)
        
        for line in infile:
            
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                column_names = line[1:].rstrip().split('\t')
                continue
                
            fields = line.rstrip().split('\t')
            if 'missense_variant' in fields[6]:
                #if 'NMD_transcript_variant' in fields[6]:
                #    continue
                line = {k: v for k, v in zip(column_names, fields)}
            else:
                continue
    
            if 'STRAND=-1' in fields[-1]:
                strand = '-'
            else:
                strand = '+'
    
            codon_pos = (int(line['CDS_position'])-1) % 3    
            altbase = line['Codons'].split('/')[1][codon_pos]
            
            altbase = complement[altbase] if strand == '-' else altbase
            
            if altbase != line['Allele']:
                g_err[line['Gene']] += 1
                continue
        
            dna = cds.fetch(line['Feature'])
            
            if len(dna) % 3:
                dna += 'NN' if (len(dna) % 3) == 1 else 'N'
            
            try:
                aaseq = Seq(dna, IUPAC.unambiguous_dna).translate()
            except:
                g_err[line['Gene']] += 1
                continue
            
            refaa, altaa = line['Amino_acids'].split('/')
            ppos = int(line['Protein_position']) - 1
            
            if str(aaseq[ppos]) != refaa:
                g_err[line['Gene']] += 1
                continue
            
            try:
                start, end, p = get_coord(ppos, aaseq)
            except Mismatch:
                g_err[line['Gene']] += 1
                continue
            
            try:
                reftrip = str(aaseq[start:end])
                if p == 1:
                    alttrip = reftrip[0] + altaa + reftrip[2]
                elif p == 0:
                    alttrip = altaa + reftrip[1] + reftrip[2]
                else:
                    alttrip = reftrip[0] + reftrip[1] + altaa
            except IndexError:
                g_err[line['Gene']] += 1
                continue
                
            cos = get_cosine(reftrip, alttrip)
            
            sift_str, sift_sc, polyphen_str, polyphen_sc = parse_scores(line['Extra'])
            
            out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                line['Uploaded_variation'],
                line['Location'],
                line['Allele'],
                line['Gene'],
                line['Feature'],
                strand,
                line['Protein_position'],
                codon_pos,
                reftrip,
                alttrip,
                cos,
                sift_str,
                sift_sc,
                polyphen_str,
                polyphen_sc
            ))
            
            i += 1
            g_succ[line['Gene']] += 1
            

pd.DataFrame.from_dict(g_err,orient='index',columns=['n_err']).to_csv(out_err)
pd.DataFrame.from_dict(g_succ,orient='index',columns=['n_succ']).to_csv(out_succ)


'''

filters output generated above, makes sure only one 3-gram is reported for every missense variant

this version of the script uses the polyphen and sift scores 
filtering works like this:

1) in case of multiple 3-grams for the same SNP and codon position, keep the one with the most common one(s)
2) If 1) results in ties, keep the ones with the highest "impact" (mean of polyphen and (1-sift) )
3) in case there are still ties (same SNP but multiple reading frames, i.e. different codon positions) keep the first entry (this should almost never (?) happen)

'''

fe = pd.read_csv(outfile, sep='\t')

locations = fe.Location.str.split(':',expand=True)
locations.loc[:,1] = locations.loc[:,1].astype(int)
locations.columns = ['chr','pos']

# "standardize" to codon start-positions
fe = fe.join(locations)
fe['pos_standardized'] = np.where(fe['Strand'] == '+', fe['pos'] - fe['codon_position'], fe['pos'] + fe['codon_position'])

fe = fe.dropna()

# should we really take mean values here?
fe_f = fe.groupby(['Gene','Uploaded_variation','Location','codon_position','ref','alt','cosine_similarity','pos_standardized']).agg({'polyphen_score':np.mean, 'sift_score':np.mean, 'pos':'size'})
fe_f.rename({'pos':'n'}, axis=1, inplace=True)
fe_f.reset_index(inplace=True)

# in case of multiple posibilities by gene, keep the most common ones
n_max = fe_f.groupby(['Gene','Uploaded_variation'])['n'].transform(np.max)
n_max = n_max == fe_f['n']
fe_f = fe_f[n_max]

# for variants with multiple consequences per gene, keep the one with the largest impact:
fe_f['impact'] = pd.concat([fe_f['polyphen_score'], 1. - fe_f['sift_score']], axis=1).mean(axis=1)
imp_max = fe_f.groupby(['Gene','Uploaded_variation'])['impact'].transform(np.max)
imp_max = imp_max == fe_f['impact']
fe_f = fe_f[imp_max]

# if multiple reading frames are present, keep the first
cod_min = fe_f.groupby(['Gene','Uploaded_variation'])['codon_position'].transform(np.min)
cod_min = cod_min == fe_f['codon_position']
fe_f = fe_f[cod_min]

fe_f = fe_f[~fe_f[['Gene','Uploaded_variation']].duplicated()]

fe_f.to_csv('{}_filtered.tsv'.format(args.out_prefix), sep='\t', index=False)
