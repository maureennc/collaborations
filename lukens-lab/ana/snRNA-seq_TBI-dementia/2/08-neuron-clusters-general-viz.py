#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 09:37:01 2024

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
plt.rcParams['font.size'] = 24

###############################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '5-neuron-recluster-hvg.h5ad'))
sc.pl.umap(adata, color = 'neuron_cluster')

## Prepare data
adata = adata[adata.obs['neuron_cluster'] != 'Inhibitory 10'].copy()
adata.obs['neuron_cluster'].replace({'Mixed Ex/Inhib': 'Mixed Ex/In'}, inplace=True)


## Cluster markers

## Save as csv
csv_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/neuron_cluster/cluster-markers"

cluster_markers = pd.read_csv(os.path.join(csv_dir, 'neuron-cluster-markers-wilcoxon.csv'))
cluster_markers

###############################################################################

# GENE LIST

order = ['Excitatory 18', 'Excitatory 2', 'Excitatory 3', 'Excitatory 8', 'Excitatory 11', 'Excitatory 4', 'Excitatory 16', 'Excitatory 1', 'Excitatory 6', 'Excitatory 12', 'Excitatory 7', 'Excitatory 9', 'Excitatory 14', 'Excitatory 17', 'Excitatory 5', 'Excitatory 13', 'Excitatory 10', 'Excitatory 15', 'Excitatory 19', 'Inhibitory 4', 'Inhibitory 7', 'Inhibitory 2', 'Inhibitory 3', 'Inhibitory 5', 'Inhibitory 9', 'Inhibitory 6', 'Inhibitory 8', 'Inhibitory 1', 'Mixed Ex/In']


# Sort cluster_markers by 'scores' and set 'group' as categorical with custom order
cluster_markers_sorted = cluster_markers.sort_values(by='scores', ascending=False)
cluster_markers_sorted['group'] = pd.Categorical(cluster_markers_sorted['group'], categories=order, ordered=True)

# Group by 'group' and select the top 3 entries for each group
genes = cluster_markers_sorted.groupby('group').head(3).reset_index(drop=True)

# Create a gene list, ordered by the 'group' category (which follows the custom order)
gene_list = []
for grp in order:
    gene_list.extend(genes[genes['group'] == grp]['names'].tolist())

# Plot the matrix plot with the reordered gene list
sc.pl.matrixplot(adata, var_names=gene_list, groupby='neuron_cluster', standard_scale='var', dendrogram=True)

sc.pl.matrixplot(adata, var_names=gene_list, groupby='neuron_cluster', dendrogram=True)

###############################################################################

# UMAP
plt.rcdefaults()

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 16


sc.pl.umap(adata, color = 'neuron_cluster', title = 'Neuronal clusters')
sc.pl.umap(adata, color = 'neuron_cluster', legend_loc = 'on data', size = 100, legend_fontoutline= 5, legend_fontsize = 25, title = 'Neuronal clusters')

###############################################################################


# CELL COMPOSITION ANALYSIS - NUMBER OF EXCITATORY AND INHIBITORY NEURONS PER GROUP


ipsilateral = ['A', 'B', 'C', 'E']
adata = adata[adata.obs['group'].isin(ipsilateral)].copy()

condition = {
    'A': 'Sham + AAV-GFP',
    'B': 'Sham + AAV-VEGFC',
    'C': 'TBI + AAV-GFP',
    'E': 'TBI + AAV-VEGFC'
}

adata.obs['condition'] = adata.obs['group'].map(condition)


###


# Convert the adata.obs to a DataFrame
df = pd.DataFrame(adata.obs)

# Group by 'condition' and 'cell_type' (assuming 'cell_type' contains Excitatory and Inhibitory labels)
condition_count = df.groupby(['condition', 'cell_type']).size().reset_index(name='Count')


# Plot the total number of cells per condition and neuron type (excitatory vs inhibitory)
plt.figure(figsize=(6, 5))
g = sns.barplot(
    x='condition', 
    y='Count', 
    hue='cell_type',  # Split by cell type (Excitatory/Inhibitory)
    data=condition_count, 
    palette='viridis'
)

# Customize plot appearance
plt.title('')  # Optionally add title
plt.xlabel('')  # Optionally customize x label
plt.ylabel('Number of Neurons')
plt.xticks(rotation=45, ha = 'right')
plt.tight_layout()

# Move the legend outside the plot
plt.legend(title='Neuron Type', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

plt.grid(False)
plt.show()


###############################################################################

# CELL COMPOSITION ANALYSIS - CELL COUNT


df = pd.DataFrame(adata.obs)

df['Cell type'] = df['cell_type'].astype('category')
df['Group'] = df['condition'].astype('category')
df['Cluster'] = df['neuron_cluster'].astype('category')

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

g._legend.set_bbox_to_anchor((1, .8))  # (x, y) The position of the legend's bounding box
g._legend.set_title('Group') 

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.show()

###############################################################################

# CELL COMPOSITION ANALYSIS - FREQUENCY
# Group by 'group' and 'neuron_cluster' to get the raw counts for each neuron cluster within each group

groups = ['A', 'B', 'C', 'E']

cluster_counts = adata.obs.groupby(['condition', 'neuron_cluster']).size().reset_index(name='Count')
print(cluster_counts)

# Calculate the total number of neurons in each group
total_neurons_per_group = adata.obs.groupby('condition').size().reset_index(name='Total_Neurons_group')
print(total_neurons_per_group)

# Merge the total neuron counts into the cluster counts DataFrame
cluster_counts = pd.merge(cluster_counts, total_neurons_per_group, on='condition')
cluster_counts

# Calculate the frequency of each neuron cluster within each group
cluster_counts['Frequency'] = cluster_counts['Count'] / cluster_counts['Total_Neurons_group']
cluster_counts
# Plot the cluster frequencies using seaborn
plt.figure(figsize=(10, 6))
g = sns.catplot(
    x='neuron_cluster', 
    y='Frequency',  # Use the calculated frequency
    hue='condition', 
    data=cluster_counts, 
    kind='bar', 
    height=5, 
    aspect=2, 
    palette='viridis'
)

# Customize the plot
g.set_axis_labels("Neuron Cluster", "% Neuronal Population per Group")  # Custom x and y labels
g._legend.set_title('Group')
g._legend.set_bbox_to_anchor((0.9, .9))  # Adjust the legend position
plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.show()


###############################################################################

# DENDROGRAM

sc.pl.dendrogram(adata, groupby = 'neuron_cluster')
