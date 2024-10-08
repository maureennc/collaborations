#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 21:55:31 2024

@author: maureen

Used in 08-neuron-cluster-de
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

## Prepare data
adata = adata[adata.obs['neuron_cluster'] != 'Inhibitory 10'].copy()
adata.obs['neuron_cluster'].replace({'Mixed Ex/Inhib': 'Mixed Ex/In'}, inplace=True)

## Extract HVG gene list
sc.pp.calculate_qc_metrics(adata, log1p = False, inplace = True)
hvg_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/neuron_cluster/cluster-markers"
hvg_list = pd.DataFrame(adata.var)
#hvg_list.to_csv(os.path.join(hvg_dir, 'hvg-list.csv'))

###############################################################################

# MARKER GENES FOR EACH CLUSTER

sc.tl.rank_genes_groups(adata, groupby='neuron_cluster', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, fontsize = 12)

cluster_markers = sc.get.rank_genes_groups_df(adata, group = None)
cluster_markers = cluster_markers[cluster_markers['pvals_adj'] < 0.05].copy()
cluster_markers

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




## Excitatory 6
exc6_genes = cluster_markers[cluster_markers['group'] == 'Excitatory 6'].copy()
exc6_genes_sorted = exc6_genes.sort_values(by='scores', ascending=False)

top_25_genes = exc6_genes_sorted['names'].head(25)
bottom_25_genes = exc6_genes_sorted['names'].tail(25).iloc[::-1]

genes = top_25_genes.tolist() + bottom_25_genes.tolist()

## Create grouping for visualization
adata.obs['Ex6_vs_rest'] = adata.obs['neuron_cluster'].apply(lambda x: 'Exc6' if x == 'Excitatory 6' else 'Rest')

sc.pl.heatmap(adata, var_names = genes, groupby = 'Ex6_vs_rest', standard_scale = 'var', use_raw = True)

sc.pl.matrixplot(adata, var_names = genes, groupby = 'Ex6_vs_rest', standard_scale = 'var', use_raw = True)


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

subset_df = sc.get.rank_genes_groups_df(subset, group = None)

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
