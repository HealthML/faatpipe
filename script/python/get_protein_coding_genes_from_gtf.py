#!/usr/bin/env python
# coding: utf-8


import pyranges
import os
import numpy as np
import logging

# gtf file with gene annotations
gtf_path = snakemake.input.gtf
target_regions_path = snakemake.input.target_regions

logging.basicConfig(filename=snakemake.log[0])


anno = pyranges.read_gtf(gtf_path)

# filter 
protein_coding_genes = anno[(anno.Feature == 'gene') & (anno.gene_biotype == 'protein_coding')]

target_regions = pyranges.read_bed(target_regions_path)
protein_coding_genes = protein_coding_genes.overlap(target_regions)

logging.info('found {} protein coding genes.'.format(len(protein_coding_genes)))

id_name = np.array(['_'.join([i, n]) for i, n in zip(protein_coding_genes.gene_id, protein_coding_genes.gene_name)])

protein_coding_genes = protein_coding_genes.drop()
protein_coding_genes.Name = id_name


os.makedirs(os.path.dirname(snakemake.output.pc_genes_bed), exist_ok=True)

protein_coding_genes.to_bed(snakemake.output.pc_genes_bed)
protein_coding_genes[protein_coding_genes.Strand == '-'].to_bed(snakemake.output.pc_genes_bed_minus)
protein_coding_genes[protein_coding_genes.Strand == '+'].to_bed(snakemake.output.pc_genes_bed_plus)

