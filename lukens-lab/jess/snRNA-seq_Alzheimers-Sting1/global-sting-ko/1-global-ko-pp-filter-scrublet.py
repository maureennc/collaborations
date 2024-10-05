#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 10:13:31 2024

@author: maureen
"""

import os
from pathlib import Path
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.rcParams['font.family'] = 'Arial'

import seaborn as sns
import pandas as pd
import numpy as np
import scanpy as sc
#sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)

from scipy.sparse import csr_matrix
import scipy
import scrublet as scr

# environment = 'sc-pp'

################################################################################################################################

# IMPORT DATA

## WT vs. Global KO
data_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/1_5xSting_Global-KO/"

adata_wt = sc.read_10x_mtx(os.path.join(data_dir, "Sample1_5xSting_globalWT/1/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
adata_ko = sc.read_10x_mtx(os.path.join(data_dir, "Sample2_5xSting_globalKO/2/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)

## Check for Sting1
print("Gene present in adata.var:", 'Tmem173' in adata_wt.var.index)
print("Gene present in adata.var:", 'Tmem173' in adata_ko.var.index)

################################################################################################################################

## ANNOTATE MITOCHONDRIAL AND RIBOSOMAL GENES

adata_list = [adata_wt, adata_ko]

for i, adata in enumerate(adata_list):
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['ribosomal'] = adata.var_names.str.match('^(Rpl|Rps)\\d+')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribosomal'], percent_top=None, log1p=False, inplace=True)


################################################################################################################################

# VISUALIZE QC DATA

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribosomal'])
    

################################################################################################################################

# FILTER CELLS BASED ON QC THRESHOLDS

## Total counts
### Floor threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before total_counts floor threshold: {adata.shape}")
    adata_list[i] = adata[adata.obs['total_counts'] > 1000, :].copy()
    print(f"Dataset {i} shape after total_counts floor threshold: {adata_list[i].shape}")


### Ceiling threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before total_counts ceiling threshold: {adata.shape}")
    total_counts_ceiling = adata.obs['total_counts'].quantile(0.95)
    adata_list[i] = adata[adata.obs['total_counts'] < total_counts_ceiling, :]
    print(f"Dataset {i} shape after total_counts ceiling threshold: {adata_list[i].shape}")

## Number of genes per cell
### Floor threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before n_genes_by_counts floor threshold: {adata.shape}")
    adata_list[i] = adata[adata.obs['n_genes_by_counts'] > 500, :].copy()
    print(f"Dataset {i} shape after n_genes_by_counts floor threshold: {adata_list[i].shape}")

### Ceiling threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before n_genes_by_counts ceiling threshold: {adata.shape}")
    ngenes_ceiling = adata.obs['n_genes_by_counts'].quantile(0.97)
    adata_list[i] = adata[adata.obs['n_genes_by_counts'] < ngenes_ceiling, :]
    print(f"Dataset {i} shape after n_genes_by_counts ceiling threshold: {adata_list[i].shape}")

## Ribosomal counts
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before pct_counts_ribosomal ceiling threshold: {adata.shape}")
    ribo_ceiling = adata.obs['pct_counts_ribosomal'].quantile(0.99)
    adata_list[i] = adata[adata.obs['pct_counts_ribosomal'] < ribo_ceiling, :]
    print(f"Dataset {i} shape after pct_counts_ribosomal ceiling threshold: {adata_list[i].shape}")

adata_wt, adata_ko = adata_list

## Mitochondrial percentage filter
for i, adata in enumerate(adata_list):
    adata_list[i] = adata[adata.obs.pct_counts_mt < 5, :].copy()

adata_wt, adata_ko = adata_list

################################################################################################################################

# RECALCULATE QC METRICS

for i, adata in enumerate(adata_list):
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribosomal'], percent_top=None, log1p=False, inplace=True)

## Print QC summary data
qc_metrics = ['total_counts', 'pct_counts_mt', 'pct_counts_ribosomal', 'n_genes_by_counts']

for i, adata in enumerate(adata_list, start=1):
    print(f"Dataset {i}:")
    for metric in qc_metrics:
        mean_value = adata.obs[metric].mean()
        print(f"Mean {metric}: {mean_value}")
    print("-" * 30) 

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribosomal'])


################################################################################################################################

# RUN SCRUBLET FOR INITIAL DOUBLET SCORING

def run_scrublet(adata_list):

    scrublet_rows = []

    for i, adata in enumerate(adata_list):
        # Run Scrublet
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.10, random_state = 0)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                                  min_cells=3, 
                                                                  min_gene_variability_pctl=85, 
                                                                  n_prin_comps=50)
        adata.obs['doublet_scores'] = doublet_scores
        adata.obs['predicted_doublets'] = predicted_doublets

        # Plot Histograms
        plt.figure(figsize=(10, 6))
        sns.histplot(scrub.doublet_scores_obs_, bins=30, color="blue", label="Observed", kde=True)
        sns.histplot(scrub.doublet_scores_sim_, bins=30, color="red", label="Simulated", kde=True)
        plt.title(f'Scrublet Doublet Score Distribution for Sample {i+1}')
        plt.xlabel('Doublet Score')
        plt.ylabel('Density')
        plt.legend()
        plt.grid(False)
        plt.show()

        # Extract cell barcodes
        cell_barcodes = adata.obs.index

        # Store Scrublet data with cell barcodes for each AnnData object in the list
        for barcode, obs_score, sim_score, pred_doublet in zip(cell_barcodes, scrub.doublet_scores_obs_, scrub.doublet_scores_sim_, predicted_doublets):
            scrublet_rows.append({'Sample_Index': i+1, 
                                  'Cell_Barcode': barcode,
                                  'Observed_Score': obs_score, 
                                  'Simulated_Score': sim_score, 
                                  'Predicted_Doublet': pred_doublet})

    # Create a DataFrame from the list of rows
    scrublet_df = pd.DataFrame(scrublet_rows)
    return scrublet_df

## Run Scrublet on list
scrublet_result = run_scrublet(adata_list)


################################################################################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/global-sting-ko/h5ad"

adata_wt.X = csr_matrix(adata_wt.X)
adata_ko.X = csr_matrix(adata_ko.X)

adata_wt.write_h5ad(os.path.join(save_dir, 'global-sting-wt-filtered.h5ad'))
adata_ko.write_h5ad(os.path.join(save_dir, 'global-sting-ko-filtered.h5ad'))

################################################################################################################################