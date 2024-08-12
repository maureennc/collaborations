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
import scrublet as scr
import scanpy as sc
import scanpy.external as sce
sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)
from scipy.sparse import csr_matrix
import scipy
from scipy.sparse import csr_matrix
#import scvi
# environment = 'sc'


import matplotlib
#matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt


os.chdir("/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq")

# Global KO

## Both files are in the KO directory
adata1 = sc.read_10x_mtx("1_5xSting_Global-KO/Sample2_5xSting_globalKO/1/filtered_feature_bc_matrix", var_names='gene_symbols', make_unique=True)
adata2 = sc.read_10x_mtx("1_5xSting_Global-KO/Sample2_5xSting_globalKO/2/filtered_feature_bc_matrix", var_names='gene_symbols', make_unique=True)



################################################################################################################################

# BASIC FILTERING

sc.pp.filter_genes(adata1, min_cells=3)
sc.pp.filter_cells(adata1, min_genes=200)
adata1.var['mt'] = adata1.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata1, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pp.filter_genes(adata2, min_cells=3)
sc.pp.filter_cells(adata2, min_genes=200)
adata2.var['mt'] = adata2.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata2, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

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


adata_list = [adata1, adata2]
scrublet_result = run_scrublet(adata_list)

def filter_doublets(adata_list):
    for adata in adata_list:
        adata._inplace_subset_obs(adata.obs['doublet_scores'] < 0.5)

filter_doublets(adata_list)

################################################################################################################################

# QC DATA VISUALIZATION

qc_metrics = ['total_counts', 'n_genes', 'n_genes_by_counts', 'total_counts_mt', 'pct_counts_mt' ]

custom_titles = {
    'total_counts': 'Total Counts Distribution',
    'n_genes': 'Number of Genes',
    'n_genes_by_counts': 'Genes by Counts',
    'total_counts_mt': 'Total Mitochondrial Counts',
    'pct_counts_mt': 'Mitochondrial Counts Percentage'
}

# Loop for QC visualizations

plt.figure(figsize=(15, 12))
for i, metric in enumerate(qc_metrics, 1): #Matplotlib's subplot indexing starts at 1, not 0
    plt.subplot(2, 3, i)
    
    # Prepare the data
    data1 = adata1.obs[metric]
    data2 = adata2.obs[metric]

    combined_df = pd.concat([
        pd.DataFrame({metric: data1, 'Dataset': 'WT'}),
        pd.DataFrame({metric: data2, 'Dataset': 'KO'})
    ])

    sns.violinplot(x='Dataset', y=metric, data=combined_df)
    sns.stripplot(x='Dataset', y=metric, data=combined_df, color='k', size = 1)
    
    ax = plt.gca()
    ax.yaxis.set_major_locator(ticker.AutoLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

    plt.title(custom_titles[metric])

plt.suptitle('QC: Global STING KO', y = 0.95) # decreasing y value moves title down
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

################################################################################################################################

# MITOCHONDRIAL FILTER

adata1 = adata1[adata1.obs.pct_counts_mt < 5, :].copy()
adata2 = adata2[adata2.obs.pct_counts_mt < 5, :].copy()

################################################################################################################################

# PRE-PROCESSING

## Normalization and Log1p transformation
adata1.layers['counts'] = adata1.X.copy()
sc.pp.normalize_total(adata1, target_sum = 1e4)
adata1.layers['normalized'] = adata1.X.copy()
sc.pp.log1p(adata1)
adata1.layers['log1p'] = adata1.X.copy()
adata1.raw = adata1


adata2.layers['counts'] = adata2.X.copy()
sc.pp.normalize_total(adata2, target_sum = 1e4)
adata2.layers['normalized'] = adata2.X.copy()
sc.pp.log1p(adata2)
adata2.layers['log1p'] = adata2.X.copy()
adata2.raw = adata2


########


## Set up and train model
scvi.model.SCVI.setup_anndata(
    adata1,
    layer='counts',
    #categorical_covariate_keys=[''],
    continuous_covariate_keys=['total_counts'])

model = scvi.model.SCVI(adata1)
model.train()

## Save model
save_dir = ('./scVI')
os.makedirs(save_dir, exist_ok=True)
model_dir = os.path.join(save_dir, "scvi_model_2024-01-11_cosmx")
model.save(model_dir, overwrite=True)

## Extract model inputs
SCVI_LATENT_KEY = "X_scVI"

latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
latent.shape

#####

### You can get back an AnnData of the object in .raw by calling .raw.to_adata()

## Regress out effects of confounding covariates
sc.pp.regress_out(adata1, ['total_counts'])
sc.pp.regress_out(adata2, ['total_counts'])

## Select HVGs
sc.pp.highly_variable_genes(adata1, n_top_genes=3000, subset=True, layer='counts', flavor="seurat_v3")
sc.pp.highly_variable_genes(adata2, n_top_genes=3000, subset=True, layer='counts', flavor="seurat_v3")

## Scale datasets
sc.pp.scale(adata1, max_value=10)
sc.pp.scale(adata2, max_value=10)


################################################################################################################################

# PCA & CLUSTERING

sc.tl.pca(adata1, n_comps = 30)
sc.pl.pca_variance_ratio(adata1, log=True, n_pcs=30)
sc.pp.neighbors(adata1, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata1)
sc.tl.leiden(adata1, resolution=0.2)


sc.tl.pca(adata2, n_comps = 30)
sc.pl.pca_variance_ratio(adata2, log=True, n_pcs=30)
sc.pp.neighbors(adata2, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata2)
sc.tl.leiden(adata2, resolution=0.2)


################################################################################################################################

# CONTROL - MARKERS

sc.tl.rank_genes_groups(adata1, groupby='leiden', use_raw = True, method = 'wilcoxon')
sc.pl.rank_genes_groups(adata1, n_genes=30)

sc.pl.umap(adata1, color = ['Rbfox3', 'Gad1', 'Slc17a6', 'Snap25',
                            'Mbp', 'Plp1', 'Sox10', 'Pdgfra', 'Ptprz1',
                            'Slc1a2', 'Apoe', 'Atp1a2', 'Aqp4', 'Aldh1l1', 'Elmo1',
                            'Hexb', 'Tgfbr1', 'C1qa',
                            'Cldn5', 'Flt1', 'Vtn',
                            'mt-Atp6', 'mt-Co1', 'Rpl13', 'Ptgds', 'Hsp90ab1', 
                            'Pecam1', 'Rgs5'])


markers1 = sc.get.rank_genes_groups_df(adata1, group = None)
markers1 = markers1[(markers1.pvals_adj < 0.05) & (markers1.logfoldchanges > .5)]
markers1


markers1.loc[markers1['group'] == '0'].head(20)
markers1.loc[markers1['group'] == '1'].head(20)
markers1.loc[markers1['group'] == '2'].head(20)
markers1.loc[markers1['group'] == '3'].head(20)
markers1.loc[markers1['group'] == '4'].head(20)
markers1.loc[markers1['group'] == '5'].head(20)
markers1.loc[markers1['group'] == '6'].head(30)
markers1.loc[markers1['group'] == '7'].head(20)
markers1.loc[markers1['group'] == '8'].head(20)
markers1.loc[markers1['group'] == '9'].head(20)
markers1.loc[markers1['group'] == '10'].head(30)
markers1.loc[markers1['group'] == '11'].head(20)
markers1.loc[markers1['group'] == '12'].head(20)
markers1.loc[markers1['group'] == '13'].head(20)
markers1.loc[markers1['group'] == '14'].head(20)


Cell_class = { 
"0": "Neuron",
"1": "Neuron", # Excitatory
"2": "Oligodendrocyte", 
"3": "Astrocyte", #Apoe+
"4": "Macrophage",  # mt-Co1, Apoe, Ptgds, Hsp90ab1, Cst3, Rpl13
"5": "Neuron", 
"6": "Ptgds+", # 
"7": "OPC",
"8": "Neuron", 
"9": "Vascular",
"10": "Neuron",
"11": "Neuron", 
"12": "Neuron",
"13": "Neuron",
"14": "Nnat+"
}

adata1.obs['Cell_class'] = adata1.obs.leiden.map(Cell_class)


sc.pl.umap(adata1, color = ['Cell_class'], legend_loc = 'on data', legend_fontsize = 10)
sc.pl.umap(adata1, color = ['leiden'], legend_loc = 'on data', legend_fontsize = 20)



################################################################################################################################

# KO - MARKERS

sc.tl.rank_genes_groups(adata2, groupby='leiden', use_raw = True, method = 'wilcoxon')
sc.pl.rank_genes_groups(adata2, n_genes=30)

sc.pl.umap(adata2, color = ['Rbfox3', 'Gad1', 'Slc17a6', 'Snap25',
                            'Mbp', 'Plp1', 'Sox10', 'Pdgfra', 'Ptprz1',
                            'Slc1a2', 'Apoe', 'Atp1a2', 'Aqp4', 'Aldh1l1', 'Elmo1',
                            'Hexb', 'Tgfbr1', 'C1qa',
                            'Cldn5', 'Flt1', 'Vtn',
                            'mt-Atp6', 'mt-Co1', 'Rpl13', 'Ptgds', 'Hsp90ab1', 
                            'Pecam1', 'Rgs5'])

sc.pl.umap(adata2, color = 'leiden')

markers2 = sc.get.rank_genes_groups_df(adata2, group = None)
markers2 = markers2[(markers2.pvals_adj < 0.05) & (markers2.logfoldchanges > .5)]
markers2


markers2.loc[markers2['group'] == '0'].head(20)
markers2.loc[markers2['group'] == '1'].head(20)
markers2.loc[markers2['group'] == '2'].head(20)
markers2.loc[markers2['group'] == '3'].head(20)
markers2.loc[markers2['group'] == '4'].head(20)
markers2.loc[markers2['group'] == '5'].head(20)
markers2.loc[markers2['group'] == '6'].head(30)
markers2.loc[markers2['group'] == '7'].head(20)
markers2.loc[markers2['group'] == '8'].head(20)
markers2.loc[markers2['group'] == '9'].head(20)
markers2.loc[markers2['group'] == '10'].head(30)
markers2.loc[markers2['group'] == '11'].head(20)
markers2.loc[markers2['group'] == '12'].head(20)
markers2.loc[markers2['group'] == '13'].head(20)
markers2.loc[markers2['group'] == '14'].head(20)
markers2.loc[markers2['group'] == '15'].head(20)
markers2.loc[markers2['group'] == '16'].head(20)
markers2.loc[markers2['group'] == '17'].head(20)
markers2.loc[markers2['group'] == '18'].head(20)
markers2.loc[markers2['group'] == '19'].head(20)
markers2.loc[markers2['group'] == '20'].head(20)
markers2.loc[markers2['group'] == '21'].head(20)



Cell_class = { 
"0": "Mito-high",
"1": "Oligodendrocytes",
"2": "Neuron", 
"3": "Astrocyte", #Apoe+
"4": "Neuron",
"5": "Neuron", 
"6": "Neuron", # 
"7": "",
"8": "", 
"9": "Macrophage",
"10": "",
"11": "", 
"12": "Neuron",
"13": "",
"14": "",
"15": "Vascular",
"16": "",
"17": "Neuron",
"18": "Neuron",
"19": "Astrocyte",
"20": "Vascular",
"21": "Neuron"
}

adata1.obs['Cell_class'] = adata1.obs.leiden.map(Cell_class)


sc.pl.umap(adata2, color = ['Cell_class'], legend_loc = 'on data', legend_fontsize = 10)
sc.pl.umap(adata2, color = ['leiden'], legend_loc = 'on data', legend_fontsize = 20)


################################################################################################################################

# Neurons
sc.pl.umap(adata1, color = ['Rbfox3', 'Snap25', 'Gad1', 'Gad2', 'Slc17a6', 'Slc17a7'], ncols = 3)
sc.pl.umap(adata2, color = ['Rbfox3', 'Snap25', 'Gad1', 'Gad2', 'Slc17a6', 'Slc17a7'], ncols = 3)

# Macrophages
sc.pl.umap(adata1, color = ['Hexb', 'P2ry12', 'Tgfbr1', 'Itgam', 'Itgax', 'Mpeg1', 'Cx3cr1', 'Mrc1', 'H2-Ab1'], ncols = 3)
sc.pl.umap(adata2, color = ['Hexb', 'P2ry12', 'Tgfbr1', 'Itgam', 'Itgax', 'Mpeg1', 'Cx3cr1', 'Mrc1', 'H2-Ab1'], ncols = 3)

# Astrocytes
sc.pl.umap(adata1, color = ['Gfap', 'Aldh1l1', 'Aqp4', 'Slc1a2', 'Atp1a2', 'Apoe'], ncols = 3)
sc.pl.umap(adata2, color = ['Gfap', 'Aldh1l1', 'Aqp4', 'Slc1a2', 'Atp1a2', 'Apoe'], ncols = 3)

# Oligodendrocytes
sc.pl.umap(adata1, color = ['Olig1','Olig2', 'Sox10', 'Pdgfra', 'Mbp', 'Mog'], ncols = 3)
sc.pl.umap(adata2, color = ['Olig1','Olig2', 'Sox10', 'Pdgfra', 'Mbp', 'Mog'], ncols = 3)

# Vascular
sc.pl.umap(adata1, color = ['Cldn5','Cdh5', 'Flt1', 'Pecam1', 'Icam1', 'Rgs5'], ncols = 3)
sc.pl.umap(adata2, color = ['Cldn5','Cdh5', 'Flt1', 'Pecam1', 'Icam1', 'Rgs5'], ncols = 3)


# Cell stress genes
adata1.var.mt[adata1.var.mt].index
adata2.var.mt[gdata.var.mt].index

sc.pl.umap(adata1, color = ['mt-Co3', 'mt-Atp6', 'mt-Co2', 'mt-Nd2', 'Rpl13', 'Rplp1', 'Ptgds', 'Hsp90ab1', 'Sod1'], ncols = 3)
sc.pl.umap(adata2, color = ['mt-Co3', 'mt-Atp6', 'mt-Co2', 'mt-Nd2', 'Rpl13', 'Rplp1', 'Ptgds', 'Hsp90ab1', 'Sod1'], ncols = 3)

################################################################################################################################



## Export annotations

dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/mc-analysis/h5ad"
os.chdir(dir)
os.getcwd()

sc.write('sting-globalWT.h5ad', adata1)
sc.write('sting-globalKO.h5ad', adata2)


################################################################################################################################



################################################################################################################################

# CONCATENATION

adata = adata1.concatenate(adata2, batch_key='genotype', batch_categories=['WT', 'KO'] )
adata.obs # concatenation adds suffix

################################################################################################################################

# PCA AND CLUSTERING

sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

## Visualize results
sc.pl.umap(adata, color = ['genotype', 'leiden'])

################################################################################################################################

# RUN HARMONY 

#sce.pp.harmony_integrate(adata, 'genotype')
#adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
#sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
#sc.tl.umap(adata)
#sc.tl.leiden(adata, resolution=0.5)

## Visualize results
#sc.pl.umap(adata, color = ['genotype', 'leiden'])

################################################################################################################################