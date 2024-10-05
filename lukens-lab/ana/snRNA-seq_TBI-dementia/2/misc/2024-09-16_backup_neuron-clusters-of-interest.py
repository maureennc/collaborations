#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 21:55:31 2024

@author: maureen
"""

import os
import pandas as pd
import scanpy as sc
import scvi
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['figure.dpi'] = 300
plt.rcParams['figure.figsize'] =(10, 5)

plt.rcParams['font.family'] = 'Arial'

###############################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '5-neuron-recluster-hvg.h5ad'))
sc.pl.umap(adata, color = 'neuron_cluster')

## Import scVI model
model_dir = os.path.join(scvi_dir, 'model_1') # HVGs
print(model_dir)
model = scvi.model.SCVI.load(model_dir, adata=adata)
model

## Add scvi latent key to obsm
SCVI_LATENT_KEY = "X_scVI"

## For model 1 (neuronal HVGs)
latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
print(f"Model latent shape: {latent.shape}")


## Prepare data
adata = adata[adata.obs['neuron_cluster'] != 'Inhibitory 10'].copy()
adata.obs['neuron_cluster'].replace({'Mixed Ex/Inhib': 'Mixed Ex/In'}, inplace=True)

## Extract HVG gene list
sc.pp.calculate_qc_metrics(adata, log1p = False, inplace = True)
hvg_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/neuron_cluster/cluster-markers"
hvg_list = pd.DataFrame(adata.var)
hvg_list.to_csv(os.path.join(hvg_dir, 'hvg-list.csv'))

###############################################################################

# MARKER GENES FOR EACH CLUSTER

sc.tl.rank_genes_groups(adata, groupby='neuron_cluster', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, fontsize = 12)

cluster_markers = sc.get.rank_genes_groups_df(adata, group = None)

## Save as csv
save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/neuron_cluster/cluster-markers"

#cluster_markers.to_csv(os.path.join(save_dir, 'neuron-cluster-markers-wilcoxon.csv'), index = False)

# Visualize for heatmap
top3_genes = cluster_markers.groupby('group').head(3)
genes = top3_genes['names'].to_list()
genes

sc.pl.heatmap(adata, var_names=genes, groupby='neuron_cluster', cmap='viridis', standard_scale='var', dendrogram = True)

###############################################################################

# CLUSTER VS REST HEATMAP

## Mixed Ex/In
mixed_genes = cluster_markers[cluster_markers['group'] == 'Mixed Ex/In'].copy()
mixed_genes_sorted = mixed_genes.sort_values(by='scores', ascending=False)

top_25_genes = mixed_genes_sorted['names'].head(25)
bottom_25_genes = mixed_genes_sorted['names'].tail(25).iloc[::-1]

genes = top_25_genes.tolist() + bottom_25_genes.tolist()

## Create grouping for visualization
adata.obs['Mixed_vs_rest'] = adata.obs['neuron_cluster'].apply(lambda x: 'Mixed' if x == 'Mixed Ex/In' else 'Rest')

sc.pl.heatmap(adata, var_names = genes, groupby = 'Mixed_vs_rest', standard_scale = 'var', use_raw = True)

sc.pl.matrixplot(adata, var_names = genes, groupby = 'Mixed_vs_rest', standard_scale = 'var', use_raw = True)




## Excitatory 10
exc10_genes = cluster_markers[cluster_markers['group'] == 'Excitatory 10'].copy()
exc10_genes_sorted = exc10_genes.sort_values(by='scores', ascending=False)

top_25_genes = exc10_genes_sorted['names'].head(25)
bottom_25_genes = exc10_genes_sorted['names'].tail(25).iloc[::-1]

genes = top_25_genes.tolist() + bottom_25_genes.tolist()

## Create grouping for visualization
adata.obs['Ex10_vs_rest'] = adata.obs['neuron_cluster'].apply(lambda x: 'Exc10' if x == 'Excitatory 10' else 'Rest')

sc.pl.heatmap(adata, var_names = genes, groupby = 'Ex10_vs_rest', standard_scale = 'var', use_raw = True)

sc.pl.matrixplot(adata, var_names = genes, groupby = 'Ex10_vs_rest', standard_scale = 'var', use_raw = True)


###############################################################################

# MIXED EXC/INHIB CLUSTER

subset = adata[adata.obs['neuron_cluster'] == 'Mixed Ex/In']
subset.X = subset.layers['counts']

## Annotate excitatory vs. Inhibitory
subset.obs['neuron_type'] = pd.Categorical(['Unknown'] * subset.shape[0], categories=['Unknown', 'Excitatory', 'Inhibitory'])

## Excitatory
exc_expr = (subset[:, 'Slc17a6'].layers['counts'] > 0).toarray().flatten()
subset.obs.loc[exc_expr, 'neuron_type'] = 'Excitatory'

## Inhibitory
inh_expr = ((subset[:, 'Gad1'].layers['counts'] > 0).toarray().flatten() | 
            (subset[:, 'Gad2'].layers['counts'] > 0).toarray().flatten())
subset.obs.loc[inh_expr, 'neuron_type'] = 'Inhibitory'


sc.pl.umap(subset, color = 'neuron_type', groups = ['Excitatory', 'Inhibitory'])


sc.tl.rank_genes_groups(subset, groupby='neuron_type', method='wilcoxon')
sc.pl.rank_genes_groups(subset, n_genes=20, sharey=False, fontsize = 20)


###############################################################################

# MIXED EXC/INHIB CLUSTER VISUALIZE PIE CHART

filtered_subset = subset[subset.obs['neuron_type'] != 'Unknown']

groups = filtered_subset.obs['group'].unique()

## Pie charts
for group in groups:
    # Subset data for the specific group
    group_data = filtered_subset[filtered_subset.obs['group'] == group]
    
    neuron_counts = group_data.obs['neuron_type'].value_counts()
    
    plt.figure(figsize=(3, 3))
    plt.pie(neuron_counts, labels=neuron_counts.index, autopct='%1.1f%%', colors=['#4CAF50', '#FF5722'])
    plt.title(f'{group}')
    plt.show()


###############################################################################

## MIXED EXC/INHIB SCVI DE

adata.obs.loc[subset.obs.index, 'neuron_type'] = subset.obs['neuron_type']
sc.pl.umap(adata, color='neuron_type', groups=['Excitatory', 'Inhibitory'])

idx1 = adata.obs['neuron_type'] == 'Excitatory'
idx2 = adata.obs['neuron_type'] == 'Inhibitory'
result = model.differential_expression(idx1=idx1, idx2=idx2)

###############################################################################

# CLUSTER VS. REST


###############################################################################


###############################################################################
###############################################################################

## Subset adata to ipsilateral groups

## Ipsi-contra
side_mapping = {
    'A': 'Ipsilateral',
    'B': 'Ipsilateral',
    'C': 'Ipsilateral',
    'D': 'Contralateral',
    'E': 'Ipsilateral',
    'F': 'Contralateral'
}

adata.obs['side'] = adata.obs['group'].map(side_mapping)
adata = adata[adata.obs['side'] == 'Ipsilateral'].copy()
#adata = adata[adata.obs['neuron_cluster'] == 'Inhibitory 2'].copy()

###############################################################################

# Characterize Mixed population (old)

## inhib2 vs rest

# Assuming you already have your DE results for 'Inhibitory 2'
de_res = sc.get.rank_genes_groups_df(adata, group=None)
de_res = de_res[de_res['group'] == 'Inhibitory 2'].copy()

# Sort by log fold change (scores) for both positive and negative values
most_positive = de_res.sort_values(by='scores', ascending=False).head(10)
most_negative = de_res.sort_values(by='scores', ascending=True).head(10)

# Combine the top 10 most positive and most negative genes
top_genes = pd.concat([most_positive, most_negative])

# Get the list of gene names for plotting
genes_to_plot = top_genes['names'].unique().tolist()

# Plot a matrix plot with the top 10 most positive and most negative genes
sc.pl.heatmap(adata, var_names=genes_to_plot, groupby='neuron_cluster', cmap='viridis', dendrogram=True)

sc.pl.matrixplot(adata, var_names=genes_to_plot, groupby='neuron_cluster', cmap='viridis', dendrogram=True, standard_scale = 'var')

###############################################################################
markers = sc.get.rank_genes_groups_df(adata, group = None)

genes = ['Usp29', 'Grid2','Lhfpl3', 'Rmst', 'Brinp3', 'Hap1', 'Cdh18', 'Camk2d', 'Pcdh7', 'Sox2', 'Rest']
sc.pl.matrixplot(adata, var_names= genes, groupby = 'neuron_cluster', dendrogram = True, )

sc.pl.heatmap(adata, var_names= genes, groupby = 'neuron_cluster', dendrogram = True, standard_scale = 'var')




###############################################################################