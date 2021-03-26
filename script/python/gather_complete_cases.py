import pandas as pd
from util.snake import clean_str

covar = pd.read_csv(snakemake.input.covar_tsv, sep='\t')

assert covar.shape == covar.dropna().shape, 'Error: covariates are not complete! remove rows with missing values before continuing.'

pheno = pd.read_csv(snakemake.input.pheno_tsv, sep='\t')

not_nan = {}
for p in pheno.columns[1:]:
    not_nan[clean_str(p)] = pheno[[pheno.columns[0],p]].dropna()[pheno.columns[0]].values

for k, v in not_nan.items():
    with open('data/covariates/complete_cases/{}.txt'.format(k), 'x') as outfile:
        for value in v:
            outfile.write('{}\n'.format(value))
            
with open('data/covariates/complete_cases/covariates.txt', 'x') as outfile:
    for value in covar.iloc[:,0].values:
         outfile.write('{}\n'.format(value))