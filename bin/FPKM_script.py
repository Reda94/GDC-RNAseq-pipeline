#!/usr/bin/env python

import pandas as pd
import sys

Gene_Lengths = pd.read_csv(sys.argv[3]+'geneLengths.txt',sep='\t')
PC_Genes = pd.read_csv(sys.argv[3]+'EnsembleIDsPCG.txt',header=None)
Counts = pd.read_csv(sys.argv[1],names=['gene_id','Count'],sep='\t')

# Merge counts with the lengths of genes
Merged = pd.merge(Counts, Gene_Lengths, on='gene_id')

# Get subset of counts that are protein coding
Protein_Coding_Counts = Counts[Counts['gene_id'].isin(PC_Genes[0].tolist())]

# Convert counts to a float
Protein_Coding_Counts['gene_id'] = Protein_Coding_Counts['Count'].astype(float)

# Get the total protein coding counts
Num_PCG = Protein_Coding_Counts['gene_id'].sum()

# Calculate FPKM for all counts
Merged['FPKM'] = Merged['Count'] * (1000000000) / (Num_PCG * Merged['aggregate_length'])

# Calculate TPM
Merged['RPK'] = Merged['Count'] / Merged['aggregate_length']
Scale = sum(Merged['RPK']) / 1000000
Merged['TPM'] = Merged['RPK'] / Scale
Merged.to_csv(sys.argv[2]+'.FPKM.txt', columns=['gene_id','FPKM'], header=False, index=False, sep='\t', mode='a')
Merged.to_csv(sys.argv[2]+'.TPM.txt', columns=['gene_id','TPM'], header=False, index=False, sep='\t', mode='a')
