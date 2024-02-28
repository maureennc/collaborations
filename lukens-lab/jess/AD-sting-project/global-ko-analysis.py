#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 09:28:54 2024

@author: maureen
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.rcParams['font.family'] = 'Arial'
import seaborn as sns

import scanpy as sc
sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)
from scipy.sparse import csr_matrix
import scipy
import squidpy as sq
import scrublet as scr
from scipy.sparse import csr_matrix

# environment = 'sc'

os.chdir("/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq")

# Global KO

## Both files are in the KO directory
gdata1 = sc.read_10x_mtx("1_5xSting_Global-KO/Sample2_5xSting_globalKO/1/filtered_feature_bc_matrix", var_names='gene_symbols', make_unique=True)
gdata2 = sc.read_10x_mtx("1_5xSting_Global-KO/Sample2_5xSting_globalKO/2/filtered_feature_bc_matrix", var_names='gene_symbols', make_unique=True)



################################################################################################################################

# BASIC FILTERING

sc.pp.filter_genes(gdata1, min_cells=3)
sc.pp.filter_cells(gdata1, min_genes=200)
gdata1.var['mt'] = gdata1.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(gdata1, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pp.filter_genes(gdata2, min_cells=3)
sc.pp.filter_cells(gdata2, min_genes=200)
gdata2.var['mt'] = gdata2.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(gdata2, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

## Already successfully filtered by 10x Cellranger

################################################################################################################################

# DOUBLET DETECTION

def run_scrublet(adata_list):
    scrublet_rows = []

    for i, adata in enumerate(adata_list):
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.10, random_state = 0)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                                  min_cells=3, 
                                                                  min_gene_variability_pctl=85, 
                                                                  n_prin_comps=50)
        adata.obs['doublet_scores'] = doublet_scores
        adata.obs['predicted_doublets'] = predicted_doublets

        plt.figure(figsize=(10, 6))
        sns.histplot(scrub.doublet_scores_obs_, bins=30, color="blue", label="Observed", kde=True)
        sns.histplot(scrub.doublet_scores_sim_, bins=30, color="red", label="Simulated", kde=True)
        plt.title(f'Doublet Score Distribution for Sample {i+1}')
        plt.xlabel('Doublet Score')
        plt.ylabel('Density')
        plt.legend()
        plt.grid(False)
        plt.show()

        cell_barcodes = adata.obs.index

        for barcode, obs_score, sim_score, pred_doublet in zip(cell_barcodes, scrub.doublet_scores_obs_, scrub.doublet_scores_sim_, predicted_doublets):
            scrublet_rows.append({'Sample_Index': i+1, 
                                  'Cell_Barcode': barcode,
                                  'Observed_Score': obs_score, 
                                  'Simulated_Score': sim_score, 
                                  'Predicted_Doublet': pred_doublet})

    scrublet_df = pd.DataFrame(scrublet_rows)
    return scrublet_df


adata_list = [gdata1, gdata2]
scrublet_result = run_scrublet(adata_list)

def filter_doublets(adata_list):
    for adata in adata_list:
        adata._inplace_subset_obs(adata.obs['doublet_scores'] < 0.5)

filter_doublets(adata_list)

################################################################################################################################

# QC DATA VISUALIZATION

qc_metrics = ['total_counts', 'pct_counts_mt', 'total_counts_mt', 'n_genes', 'n_genes_by_counts']

custom_titles = {
    'total_counts': 'Total Counts Distribution',
    'pct_counts_mt': 'Mitochondrial Counts Percentage',
    'total_counts_mt': 'Total Mitochondrial Counts',
    'n_genes': 'Number of Genes',
    'n_genes_by_counts': 'Genes by Counts'
}

# Loop for QC visualizations

plt.figure(figsize=(15, 12))
for i, metric in enumerate(qc_metrics, 1): #Matplotlib's subplot indexing starts at 1, not 0
    plt.subplot(2, 3, i)
    
    # Prepare the data
    data1 = gdata1.obs[metric]
    data2 = gdata2.obs[metric]

    combined_df = pd.concat([
        pd.DataFrame({metric: data1, 'Dataset': 'WT'}),
        pd.DataFrame({metric: data2, 'Dataset': 'KO'})
    ])

    sns.violinplot(x='Dataset', y=metric, data=combined_df)

    ax = plt.gca()
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

    plt.title(custom_titles[metric])

plt.suptitle('QC: Global STING KO')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

################################################################################################################################

# MITOCHONDRIAL FILTER

gdata1 = gdata1[gdata1.obs.pct_counts_mt < 5, :].copy()
gdata2 = gdata2[gdata2.obs.pct_counts_mt < 5, :].copy()

################################################################################################################################

# PRE-PROCESSING

## Normalization and Log1p transformation
gdata1.layers['counts'] = gdata1.X.copy()
sc.pp.normalize_total(gdata1, target_sum = 1e4)
gdata1.layers['normalized'] = gdata1.X.copy()
sc.pp.log1p(gdata1)
gdata1.layers['log1p'] = gdata1.X.copy()
gdata1.raw = gdata1


gdata2.layers['counts'] = gdata2.X.copy()
sc.pp.normalize_total(gdata2, target_sum = 1e4)
gdata2.layers['normalized'] = gdata2.X.copy()
sc.pp.log1p(gdata2)
gdata2.layers['log1p'] = gdata2.X.copy()
gdata2.raw = gdata2


### You can get back an AnnData of the object in .raw by calling .raw.to_adata()

## Regress out effects of confounding variables
sc.pp.regress_out(gdata1, ['total_counts', 'pct_counts_mt'])
sc.pp.regress_out(gdata2, ['total_counts', 'pct_counts_mt'])

## Select HVGs
sc.pp.highly_variable_genes(gdata1, n_top_genes=3000, subset=True, layer='counts', flavor="seurat_v3")
sc.pp.highly_variable_genes(gdata2, n_top_genes=3000, subset=True, layer='counts', flavor="seurat_v3")


## Scale datasets
sc.pp.scale(gdata1, max_value=10)
sc.pp.scale(gdata2, max_value=10)


################################################################################################################################

# CONCATENATION

adata = gdata1.concatenate(gdata2, batch_key='genotype', batch_categories=['WT', 'KO'] )
adata.obs # concatenation adds suffix

################################################################################################################################

################################################################################################################################