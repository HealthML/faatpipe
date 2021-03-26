

import pandas as pd
import logging

logging.basicConfig(filename=snakemake.log[0])

vepdf = pd.read_csv(snakemake.input.tsv, sep='\t', comment='#', header=None,
                    names=['Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Consequence'])

vepdf_summary = vepdf.groupby(['Uploaded_variation', 'Location', 'Allele', 'Gene']).size()
vepdf_summary = vepdf_summary.reset_index().rename(columns={0: 'n_affected'})

vepdf_summary.to_csv(snakemake.output.tsv, index=False, sep='\t')
