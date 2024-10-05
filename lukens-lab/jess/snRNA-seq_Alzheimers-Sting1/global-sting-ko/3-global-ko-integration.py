#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 11:40:21 2024

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
sc.set_figure_params(scanpy = True, dpi = 100, dpi_save = 400, fontsize = 14, figsize = None)

################################################################################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/global-sting-ko/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, 'global-sting-ko-pp-concat-hvg.h5ad'))

################################################################################################################################

# SET UP AND TRAIN MODEL A

## Model A
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    categorical_covariate_keys=['group'],
    continuous_covariate_keys=['total_counts', 'pct_counts_mt', 'pct_counts_ribosomal'],
)

model_A = scvi.model.SCVI(adata)

scvi.train.Trainer(accelerator='cpu', devices=1)
model_A.train()

################################################################################################################################

# SAVE / IMPORT MODEL

scvi_dir = '/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/global-sting-ko/scvi'

## Save 
model_A_dir = os.path.join(scvi_dir, 'model_A')
print(model_A_dir)
#model_A.save(model_A_dir)


## Import model
scvi_dir = '/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/global-sting-ko/scvi'

model_A_dir = os.path.join(scvi_dir, 'model_A')
print(model_A_dir)
model_A = scvi.model.SCVI.load(model_A_dir, adata=adata)
model_A

################################################################################################################################

# EVALUATE TRAINED MODEL A

## Extract dictionary
training_history = model_A.history
training_history


training_history_df = pd.DataFrame(index=training_history['kl_weight'].index)

for key, df in training_history.items():
    training_history_df = training_history_df.join(df, how='outer')

## Visualize results
training_history_df.reset_index(inplace=True)

plt.figure(figsize=(5, 20))
## ELBO
plt.subplot(3, 1, 1)
plt.plot(training_history_df['epoch'], training_history_df['elbo_train'], label='ELBO')
plt.xlabel('Epochs')
plt.ylabel('ELBO')
plt.title('ELBO over Training Epochs')
plt.legend()

## Training Loss
plt.subplot(3, 1, 2)
plt.plot(training_history_df['epoch'], training_history_df['train_loss_epoch'], label='Training Loss')
plt.xlabel('Epochs')
plt.ylabel('Training Loss')
plt.title('Training Loss over Epochs')
plt.legend()

## KL Divergence (Local)
plt.subplot(3, 1, 3)
plt.plot(training_history_df['epoch'], training_history_df['kl_local_train'], label='KL Divergence (Local)')
plt.xlabel('Epochs')
plt.ylabel('KL Divergence (Local)')
plt.title('KL Divergence over Epochs')
plt.legend()

## Adjust layout
plt.tight_layout()
plt.show()


################################################################################################################################

# EXTRACT LATENT REPRESENTATION

## Check for scVI entries in obsm
adata.obsm

## add scvi latent key to obsm
SCVI_LATENT_KEY = "X_scVI"

latent = model_A.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
latent.shape

## Add scvi normalized counts layer
adata.layers['scvi_normalized'] = model_A.get_normalized_expression()
adata.layers

################################################################################################################################

# INITIAL CLUSTERING

sc.pp.neighbors(adata, use_rep = 'X_scVI', random_state = 0)
sc.tl.umap(adata, min_dist = 0.3)
sc.tl.leiden(adata, key_added = 'leiden_scVI', resolution = 1)

sc.pl.umap(adata, color = ['group', 'leiden_scVI'])

################################################################################################################################

# GET CLUSTER GENE MARKERS

sc.tl.rank_genes_groups(adata, groupby='leiden_scVI', method='wilcoxon', use_raw = True)
sc.pl.rank_genes_groups(adata, n_genes = 30)
markers = sc.get.rank_genes_groups_df(adata, None)
#markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)] # Keep all results
markers[markers['group'] == '21'].sort_values(by='logfoldchanges', ascending=False).head(50)



################################################################################################################################

# CLUSTER VALIDATION


genes = ['Rbfox3', 'Map2', 'Syn1', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 
         'Cx3cr1', 'Hexb', 'Mrc1',
         'Slc1a2', 'Slc1a3', 'Gfap', 'Aldh1l1', 'Aqp4', 'Atp1a2',
         'Pecam1', 'Tie1', 'Cldn5', 'Flt1',
         'Pdgfrb', 'Rgs5', 
         'Pdgfra', 'Cspg4', 'Sox10',
         'Mbp', 'Plp1', 'Mog'] 

sc.pl.dotplot(adata, genes, groupby = 'leiden_scVI', standard_scale = 'var', use_raw = True, dendrogram = True)


################################################################################################################################

# ANNOTATIONS

cell_type= { 
"0": "Inhibitory neuron",
"1": "Oligodendrocyte",
"2": "Astrocyte",
"3": "Excitatory neuron",
"4": "Excitatory neuron",
"5": "Excitatory neuron",
"6": "Excitatory neuron",
"7": "Inhibitory neuron",
"8": "Microglia",
"9": "Inhibitory neuron",
"10": "OPC",
"11": "Excitatory neuron",
"12": "Excitatory neuron",
"13": "Excitatory neuron",
"14": "Inhibitory neuron",
"15": "Endothelial cell",
"16": "Excitatory neuron",
"17": "Inhibitory neuron",
"18": "Excitatory neuron",
"19": "Excitatory neuron",
"20": "Excitatory neuron",
"21": "Unknown",
"22": "Inhibitory neuron",
"23": "Inhibitory neuron",
"24": "Excitatory neuron",
"25": "Excitatory neuron",
"26": "Unknown",
"27": "Inhibitory neuron",
"28": "Excitatory neuron",
"29": "Inhibitory neuron",
"30": "Excitatory neuron",
"31": "Unknown",
"32": "Other neuron",
"33": "Excitatory neuron",
"34": "Unknown",
"35": "Excitatory neuron",
"36": "Excitatory neuron"
}

adata.obs['cell_type'] = adata.obs.leiden_scVI.map(cell_type)
sc.pl.dotplot(adata, genes, groupby = 'cell_type', standard_scale = 'var', use_raw = True, dendrogram = True)


cluster= { 
"0": "Inhibitory neuron 1",
"1": "Oligodendrocyte",
"2": "Astrocyte",
"3": "Excitatory neuron 1",
"4": "Excitatory neuron 2",
"5": "Excitatory neuron 3",
"6": "Excitatory neuron 4",
"7": "Inhibitory neuron 2",
"8": "Microglia",
"9": "Inhibitory neuron 3",
"10": "OPC",
"11": "Excitatory neuron 5",
"12": "Excitatory neuron 6",
"13": "Excitatory neuron 7",
"14": "Inhibitory neuron 4",
"15": "Endothelial cell",
"16": "Excitatory neuron 8",
"17": "Inhibitory neuron 5",
"18": "Excitatory neuron 9",
"19": "Excitatory neuron 10",
"20": "Excitatory neuron 11",
"21": "Unknown 1",
"22": "Inhibitory neuron 6",
"23": "Inhibitory neuron 7",
"24": "Excitatory neuron 12",
"25": "Excitatory neuron 13",
"26": "Unknown 2",
"27": "Inhibitory neuron 8",
"28": "Excitatory neuron 14",
"29": "Inhibitory neuron 9",
"30": "Excitatory neuron 15",
"31": "Unknown 3",
"32": "Other neuron",
"33": "Excitatory neuron 16",
"34": "Unknown 4",
"35": "Excitatory neuron 17",
"36": "Excitatory neuron 18"
}

adata.obs['cluster'] = adata.obs.leiden_scVI.map(cluster)

################################################################################################################################

# VISUALIZATION

sc.pl.dotplot(adata, genes, groupby = 'cluster', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(adata, genes, groupby = 'cluster', standard_scale = 'var', use_raw = True, dendrogram = True)

sc.pl.umap(adata, color = ['leiden_scVI'], legend_loc = 'on data', legend_fontsize = 10)
sc.pl.umap(adata, color = ['cell_type'])
sc.pl.umap(adata, color = ['group'])
sc.pl.umap(adata, color = ['pct_counts_mt'])
sc.pl.umap(adata, color = ['pct_counts_ribosomal'])


################################################################################################################################

# CELL COMPOSITION ANALYSIS

df = pd.DataFrame(adata.obs)

df['Cell type'] = df['cell_type'].astype('category')
df['Group'] = df['dataset'].astype('category')
df['Cluster'] = df['cluster'].astype('category')

## Cell types
cell_type_count = df.groupby(['Group', 'Cell type']).size().reset_index(name='Count')

g = sns.catplot(
    x='Cell type',
    y='Count',
    hue='Group',  
    data=cell_type_count,
    kind='bar',
    height=5,
    aspect=2,
    palette='viridis',
    legend= True
)

g._legend.set_bbox_to_anchor((1, .9))  # (x, y) The position of the legend's bounding box
g._legend.set_title('Group') 

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.show()


## Clusters
cluster_count = df.groupby(['Group', 'Cluster']).size().reset_index(name='Count')

g = sns.catplot(
    x='Cluster',
    y='Count',
    hue='Group',  
    data=cluster_count,
    kind='bar',
    height=5,
    aspect=2,
    palette='viridis',
    legend= True
)

g._legend.set_bbox_to_anchor((1, .9))  # (x, y) The position of the legend's bounding box
g._legend.set_title('Group') 

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.show()

################################################################################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/global-sting-ko/h5ad"

adata.X = csr_matrix(adata.X)

adata.write_h5ad(os.path.join(save_dir, 'global-ko-trained-annotated.h5ad'))

################################################################################################################################


umap_genes = ['leiden_scVI', 'Mag', 'Mog', 'Aqp4', 'Atp1a2', 'Gfap', 'Hexb', 'Rbfox3', 'Slc17a6',
         'Slc17a7', 'Gad1', 'Gad2', 'Cdh5', 'Rgs5', 'Tie1', 'Pdgfra']

sc.pl.umap(adata, color = umap_genes, legend_loc = 'on data', legend_fontsize = 16)

################################################################################################################################