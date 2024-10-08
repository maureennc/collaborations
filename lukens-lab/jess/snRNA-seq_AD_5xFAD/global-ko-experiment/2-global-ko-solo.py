#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 10:29:15 2024

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
from scipy.sparse import csr_matrix

import random
import torch
import scvi

print(sns.__version__)
print(pd.__version__)
print(np.__version__)
print(sc.__version__)
print(scvi.__version__)

################################################################################################################################

# SETTINGS

## Random seed
random.seed(0)
torch.manual_seed(0)
np.random.seed(0)
scvi.settings.seed = 0

## Matplotlib
%matplotlib qt5
plt.rcParams['font.family'] = 'Arial'

## Scanpy
#sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)

################################################################################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/global-sting-ko/h5ad"

## adata used for doublet training
adata_wt = sc.read_h5ad(os.path.join(data_dir, 'global-wt-filtered.h5ad'))
adata_ko = sc.read_h5ad(os.path.join(data_dir, 'global-ko-filtered.h5ad'))

## bdata used for processing and exports
bdata_wt = sc.read_h5ad(os.path.join(data_dir, 'global-wt-filtered.h5ad'))
bdata_ko = sc.read_h5ad(os.path.join(data_dir, 'global-ko-filtered.h5ad'))

################################################################################################################################

# PREPARE DATA FOR TRAINING

## Dimensionality reduction
sc.pp.highly_variable_genes(adata_wt, n_top_genes = 3000, subset = True, flavor = 'seurat_v3')
sc.pp.highly_variable_genes(adata_ko, n_top_genes = 3000, subset = True, flavor = 'seurat_v3')

## Made adata.X CSR
adata_wt.X = csr_matrix(adata_wt.X)
adata_ko.X = csr_matrix(adata_ko.X)

################################################################################################################################

# DOUBLET DETECTION - WT

## WT - Train vae model
scvi.model.SCVI.setup_anndata(adata_wt)
vae_wt = scvi.model.SCVI(adata_wt)
vae_wt.train()

## WT - Train SOLO model
solo_wt = scvi.external.SOLO.from_scvi_model(scvi_model = vae_wt,
                                             adata = adata_wt,
                                             doublet_ratio = 3)
solo_wt.train()

## WT - Extract results
solo_wt_results = solo_wt.predict()
solo_wt_results['prediction'] = solo_wt.predict(soft = False)

print(solo_wt_results.prediction.value_counts())

################################################################################################################################

# DOUBLET DETECTION - KO

## KO - Train vae model
scvi.model.SCVI.setup_anndata(adata_ko)
vae_ko = scvi.model.SCVI(adata_ko)
vae_ko.train()

## KO - Train SOLO model
solo_ko = scvi.external.SOLO.from_scvi_model(scvi_model = vae_ko,
                                             adata = adata_ko,
                                             doublet_ratio = 3)
solo_ko.train()

## KO - Extract results
solo_ko_results = solo_ko.predict()
solo_ko_results['prediction'] = solo_ko.predict(soft = False)

print(solo_ko_results.prediction.value_counts())

################################################################################################################################

# MAP DOUBLETS TO BDATA

## Add annotations from results table
adata_wt.obs[['doublet', 'singlet', 'doublet_predictions']] = solo_wt_results[['doublet', 'singlet', 'prediction']]
adata_ko.obs[['doublet', 'singlet', 'doublet_predictions']] = solo_ko_results[['doublet', 'singlet', 'prediction']]

## Define function
def map_metrics(source_adata, target_adata):    
    for column in ['doublet', 'singlet', 'doublet_predictions']:
        target_adata.obs[column] = 'unknown' # Ensure that  target adata has columns initialized

    ## Find intersection of indices to map common barcodes
    common_barcodes = target_adata.obs_names.intersection(source_adata.obs_names)
    target_adata.obs.loc[common_barcodes, ['doublet', 'singlet', 'doublet_predictions']] = source_adata.obs.loc[common_barcodes, ['doublet', 'singlet', 'doublet_predictions']]

## Apply mapping
map_metrics(adata_wt, bdata_wt)
map_metrics(adata_ko, bdata_ko)

## Check
print(bdata_wt.obs[['doublet', 'singlet', 'doublet_predictions']].head())
print(bdata_ko.obs[['doublet', 'singlet', 'doublet_predictions']].head())


################################################################################################################################

# VISUALIZE DOUBLET PREDICTIONS

## WT
sc.pp.normalize_total(adata_wt)
sc.pp.log1p(adata_wt)
sc.tl.pca(adata_wt, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_wt, log=True, n_pcs = 50)
sc.pp.neighbors(adata_wt, n_pcs = 30)
sc.tl.umap(adata_wt)
sc.tl.leiden(adata_wt, resolution = 0.5)

sc.pl.umap(adata_wt, color = ['leiden'], legend_loc = 'on data')
sc.pl.umap(adata_wt, color = ['doublet_predictions'])


## KO
sc.pp.normalize_total(adata_ko)
sc.pp.log1p(adata_ko)
sc.tl.pca(adata_ko, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_ko, log=True, n_pcs = 50)
sc.pp.neighbors(adata_ko, n_pcs = 30)
sc.tl.umap(adata_ko)
sc.tl.leiden(adata_ko, resolution = 0.5)

sc.pl.umap(adata_ko, color = ['leiden'], legend_loc = 'on data')
sc.pl.umap(adata_ko, color = ['doublet_predictions'])


################################################################################################################################

# FILTER DOUBLETS

print(bdata_wt.obs['n_genes_by_counts'].median())
print(bdata_ko.obs['n_genes_by_counts'].median())
print(bdata_wt.obs['total_counts'].median())
print(bdata_ko.obs['total_counts'].median())

bdata_wt = bdata_wt[bdata_wt.obs['doublet_predictions'] == 'singlet'].copy()
bdata_ko = bdata_ko[bdata_ko.obs['doublet_predictions'] == 'singlet'].copy()

print(bdata_wt.obs['n_genes_by_counts'].median())
print(bdata_ko.obs['n_genes_by_counts'].median())
print(bdata_wt.obs['total_counts'].median())
print(bdata_ko.obs['total_counts'].median())

################################################################################################################################

# CONCATENATION

## Annotate groups and methods
bdata_wt.obs['group'] = 'WT'
bdata_ko.obs['group'] = 'Sting-KO'


## Perform concatenation
bdata = bdata_wt.concatenate(bdata_ko, batch_key='dataset', batch_categories=['global_wt', 'global_ko'] )
bdata.obs

################################################################################################################################

# FINISH PRE-PROCESSING

## Filter genes appearing in few cells
sc.pp.filter_genes(bdata, min_cells=3)

## Finish pre-processing (Normalization, log1p-transformation, scaling)
bdata.layers['counts'] = bdata.X.copy()
sc.pp.normalize_total(bdata)
bdata.layers['normalized'] = bdata.X.copy()
sc.pp.log1p(bdata)
bdata.layers['log1p'] = bdata.X.copy()
bdata.raw = bdata
sc.pp.scale(bdata, max_value=10)
bdata.layers['scaled'] = bdata.X.copy()

################################################################################################################################

# PRELIMINARY SELECTION OF GENES OF INTEREST AND HVGs

## Import genes of interest list
file_path = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/gene-lists/genes-of-interest.csv"
gene_list = pd.read_csv(file_path)

## Initialize the 'gene_of_interest' column to False
bdata.var['gene_of_interest'] = False

## Expand gene list with wildcards and annotate bdata.var 'gene_of_interest' column
expanded_rows = []
for _, row in gene_list.iterrows():
    if pd.isna(row['expanded']):  # If there's no expansion, include the gene directly
        if row['gene'] in bdata.var_names:
            bdata.var.loc[row['gene'], 'gene_of_interest'] = True
    else:
        expansions = row['expanded'].split(',')
        for exp in expansions:
            expanded_gene = f"{row['gene'].rstrip('*')}{exp}"
            if expanded_gene in bdata.var_names:
                bdata.var.loc[expanded_gene, 'gene_of_interest'] = True

## Select Highly Variable Genes (HVGs)
sc.pp.highly_variable_genes(bdata, n_top_genes=3000, subset=False, layer='counts', flavor="seurat_v3")
bdata.var.loc[bdata.var['highly_variable'], 'gene_of_interest'] = True

## Subset bdata to include only genes of interest
bdata_subset = bdata[:, bdata.var['gene_of_interest']].copy()
print(bdata_subset.var[['highly_variable', 'gene_of_interest']].value_counts())


################################################################################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/global-sting-ko/h5ad"

bdata_subset.X = csr_matrix(bdata_subset.X)
bdata.X = csr_matrix(bdata.X)

bdata_subset.obs['doublet'] = bdata_subset.obs['doublet'].astype(float)
bdata_subset.obs['singlet'] = bdata_subset.obs['singlet'].astype(float)

bdata.obs['doublet'] = bdata.obs['doublet'].astype(float)
bdata.obs['singlet'] = bdata.obs['singlet'].astype(float)

bdata_subset.write_h5ad(os.path.join(save_dir, 'global-sting-ko-pp-concat-hvg.h5ad'))
bdata.write_h5ad(os.path.join(save_dir, 'global-sting-ko-pp-concat-full.h5ad'))

################################################################################################################################


