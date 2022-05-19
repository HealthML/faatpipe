import pandas as pd
import numpy as np
from util.snake import clean_str

covar = pd.read_csv(snakemake.input.covar_tsv, sep='\t')

assert covar.shape == covar.dropna().shape, 'Error: covariates are not complete! remove rows with missing values before continuing.'

pheno = pd.read_csv(snakemake.input.pheno_tsv, sep='\t')

phenos = [ clean_str(x) for x in pheno.columns.tolist()[1:] ]

for i, k in enumerate(phenos):
    
    iid = pheno.iloc[:,[0,i+1]].dropna().iloc[:,0].values
    iid = np.intersect1d(iid, covar.iloc[:,0].values)
    
    with open('data/covariates/complete_cases/{}.txt'.format(k), 'x') as outfile:
        for value in iid:
            outfile.write('{}\n'.format(value))
            
with open('data/covariates/complete_cases/covariates.txt', 'x') as outfile:
    for value in covar.iloc[:,0].values:
         outfile.write('{}\n'.format(value))