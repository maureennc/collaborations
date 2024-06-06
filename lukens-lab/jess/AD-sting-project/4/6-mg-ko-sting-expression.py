#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 14:57:28 2024

@author: maureen
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.stats
from scipy.sparse import csr_matrix
import random


################################################################################################################################

# SETTINGS

## Random seed
random.seed(0)
np.random.seed(0)

## Matplotlib
%matplotlib inline
plt.rcParams['font.family'] = 'Arial'

## Scanpy
sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)

################################################################################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/mg-sting-ko/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, 'mg-model_B-cleaned-annotated.h5ad'))

adata.obs['group'] = adata.obs['group'].replace('Sting-KO', 'Sting1-KO')

################################################################################################################################

# ANNOTATE CELLS WITH NON-ZERO TMEM173 TRANSCRIPT COUNT

## Extract 'Tmem173' counts from the 'counts' layer
Tmem173_counts = adata[:, 'Tmem173'].layers['counts']

if scipy.sparse.issparse(Tmem173_counts):  # Check if the data is sparse
    Tmem173_data = Tmem173_counts.toarray()  # Convert sparse matrix to dense array
else:
    Tmem173_data = Tmem173_counts  # Use as is if it's already dense

## Generate a mask for non-zero values
non_zero_Tmem173_mask = (Tmem173_data > 0).flatten()

## Ensure the mask is properly aligned with the number of observations
if len(non_zero_Tmem173_mask) != adata.n_obs:
    raise ValueError("The length of the mask does not match the number of observations in the AnnData object.")

## Update the AnnData object with new annotations
adata.obs['Tmem173_positive'] = non_zero_Tmem173_mask
print(adata.obs['Tmem173_positive'].value_counts())

# Map boolean values to labels
adata.obs['Tmem173_expression'] = adata.obs['Tmem173_positive'].map({True: 'Tmem173-positive', False: 'Tmem173-negative'})
print(adata.obs['Tmem173_expression'].value_counts())

################################################################################################################################

# VISUALIZE DISTRIBUTION OF TMEM173+ CELLS

tmem173_positive_data = adata.obs[adata.obs['Tmem173_positive']]

# Plotting
plt.figure(figsize=(4, 3), dpi=300)
ax = sns.countplot(x='cell_type', data=tmem173_positive_data, palette='viridis')
ax.set_title('Distribution of Tmem173-positive Cells')
ax.set_xlabel('')
ax.set_ylabel('# Cells')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
ylims = plt.gca().get_ylim()
plt.yticks(np.arange(ylims[0], ylims[1], step=10))

plt.grid(False)
plt.show()

################################################################################################################################

# GROUPED BAR PLOT

plt.figure(figsize=(8, 5), dpi=300)

ax = sns.countplot(
    x='cell_type', 
    hue='group',
    data=adata.obs[adata.obs['Tmem173_positive']],
    palette='viridis',
    hue_order=['WT', 'Sting1-KO']
)

ax.set_title('Cell Type Distribution of Tmem173-positive Cells by Genotype')
ax.set_xlabel('Cell Type')
ax.set_ylabel('# Cells')
plt.xticks(rotation=45, ha='right')
ax.yaxis.grid(False)
plt.tight_layout()
plt.show()

################################################################################################################################

# Matrixplot

sc.pl.matrixplot(adata, var_names = ['Tmem173'], groupby = 'cell_type')
