#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import subprocess
import os
import hashlib
import sys
import gzip
from argparse import ArgumentParser
import logging

'''
generates a MAF report (by calling plink)
'''

def get_args():
    p = ArgumentParser()
    # meaning of the different parameters explained below
    p.add_argument('-i','--iid', help='file containing the individual-IDs that we want to include in the calculation', default=None)
    p.add_argument('-b','--bed', required=True)
    p.add_argument('-o','--out_prefix', required=True)
    p.add_argument('--plink_path', required=False, default='bin/plink', help='path to plink 1.9 binary')
    p.add_argument('--plink2_path', required=False, default='bin/plink2', help='path to plink 2 binary')
    p.add_argument('--hwe', required=False, default=0.00001, type=float)
    p.add_argument('--log', required=False, default='log.txt')
    p.add_argument('--threads', default=4, type=int)
    args = p.parse_args()
    return(args)

def read_ids(path):
    if path.endswith('.gz'):
        with gzip.open(path,'rt') as f:
            ids = [ l.rstrip() for l in f ]
    else:
        with open(path, 'r') as f:
            ids = [ l.rstrip() for l in f ]
    return ids


def write_ids(ids, path):
    with open(path, 'w') as outfile:
        for i in ids:
            outfile.write('{} {}\n'.format(str(i), str(i))) 


def main():

    # arguments
    args = get_args()
    
    if args.log is not None:
        logging.basicConfig(filename=args.log)
    else:
        logging.basicConfig()
    
    
    if args.iid is not None:
        # read the individual files
        iids = np.unique(read_ids(args.iid))
            # temporary files for plink input
        tmpfilename_iid = hashlib.sha256((''.join(iids)+args.out_prefix).encode('utf-8')).hexdigest()
        tmpfilename_iid = '{}_{}_iid.txt'.format(args.out_prefix, tmpfilename_iid)
            
        write_ids(iids, tmpfilename_iid)
        
    
    logging.info(subprocess.check_output([args.plink_path, '--version']).decode('utf-8').strip())
    
    # plink command for pre-filtering based on HWE:
    if args.hwe > 0.:
        tmpfilename_vid = '{}.snplist'.format(args.out_prefix)
        if args.iid is not None:
            plink_command = [args.plink2_path, '--keep', tmpfilename_iid, '--bfile', args.bed, '--threads', str(args.threads), '--hwe', str(args.hwe), '--write-snplist', '--out', args.out_prefix]
        else:
            plink_command = [args.plink2_path, '--bfile', args.bed, '--threads', str(args.threads), '--hwe', str(args.hwe), '--write-snplist', '--out', args.out_prefix]
        logging.info('plink command: {}'.format(' '.join(plink_command)))
        _ = subprocess.run(plink_command, check=True)

    # plink command to create count reports
    plink_command = [args.plink_path, '--bfile', args.bed]
    if args.iid is not None:
        plink_command += ['--keep', tmpfilename_iid, '--freq', 'counts', 'gz', '--threads', str(args.threads), '--out', args.out_prefix]
    else:
        plink_command += ['--freq', 'counts', 'gz', '--threads', str(args.threads), '--out', args.out_prefix ]
        
    if args.hwe > 0.:
        plink_command += ['--extract', tmpfilename_vid]
        
    logging.info('plink command: {}'.format(' '.join(plink_command)))
    
    # run plink command
    _ = subprocess.run(plink_command, check=True)

    if args.iid is not None:
        os.remove(tmpfilename_iid)
    if args.hwe > 0.:
        os.remove(tmpfilename_vid)

    # print the plink log
    with open(args.out_prefix + '.log', 'r') as infile:
        lines = infile.readlines()
        logging.info(lines)
    
    outfilepath = args.out_prefix + '.frq.counts.gz'
    counts = pd.read_csv(outfilepath, delim_whitespace=True, header=0)
    
    counts.loc[:,'CHR'] = counts['CHR'].astype(str)

    # get major / minor allele counts as they might be "flipped"
    counts['Major'] = np.where(counts['C2'] < counts['C1'], counts['C1'], counts['C2'])
    counts['Minor'] = np.where(counts['C2'] < counts['C1'], counts['C2'], counts['C1'])
    
    counts['alt_greater_ref'] = (counts['C2'] < counts['C1']).astype(int)

    n_subj = (counts.iloc[0].Major + counts.iloc[0].Minor + counts.iloc[0].G0 * 2) / 2

    # set filters
    missing = (counts.Minor == 0).values
    print('{} unobserved variants.'.format(np.sum(missing)))

    counts['MAF'] = counts.Minor / ( n_subj * 2 )
    counts['MISSING'] = counts.G0 / n_subj
    
    # export MAF report
    counts.drop(['CHR','A1','A2','C1','C2'], axis=1, inplace=True)
    counts.to_csv(args.out_prefix+'.tsv.gz', sep='\t', index=False)
    
    
    
if __name__ == '__main__':
    main()