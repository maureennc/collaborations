#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 17:05:23 2024

@author: maureen
"""

import os
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from functools import reduce
from functools import reduce

plt.rcParams['figure.dpi'] = 600
plt.rcParams['font.family'] = 'Arial'

###############################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/3/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '6-tbi-annotated-full.h5ad'))

bdata = sc.read_h5ad(os.path.join('/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/h5ad/5-neuron-recluster-full.h5ad'))

###############################################################################

# Figure 4 - UMAPs

## Change directory for scanpy save figures to work
os.chdir('/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/submission_01')

## With legend
sc.pl.umap(adata, color = ['cell_type'], legend_fontsize = 14, title = 'Cell types', save = '_cell_type.pdf')
sc.pl.umap(bdata, color = ['neuron_cluster'], legend_fontsize = 10, title = 'Neuronal clusters', save = '_neuron_cluster.pdf')

## No legend
sc.pl.umap(adata, color=['cell_type'], legend_loc=None, title='Cell types', save='_cell_type_no_legend.pdf')

sc.pl.umap(bdata, color=['neuron_cluster'], legend_loc=None, title='Neuronal clusters', save='_neuron_cluster_no_legend.pdf')

###############################################################################

# NEURON CLUSTER COMPOSITION

groups = ['A', 'B', 'C', 'E']
cdata = bdata[bdata.obs['group'].isin(groups)].copy()

## Rename group annotations

### Italicized
#condition = {
#    'A': r'Sham + AAV$^{GFP}$',
#    'B': r'Sham + AAV$^{VEGFC}$',
#    'C': r'TBI + AAV$^{GFP}$',
#    'E': r'TBI + AAV$^{VEGFC}$'
#}


### Not italicized
condition = {
    'A': r'Sham + AAV$^{\mathrm{GFP}}$',
    'B': r'Sham + AAV$^{\mathrm{VEGFC}}$',
    'C': r'TBI + AAV$^{\mathrm{GFP}}$',
    'E': r'TBI + AAV$^{\mathrm{VEGFC}}$'
}

cdata.obs['condition'] = cdata.obs['group'].map(condition)

   
## Bar chart  

plt.rcParams['font.size'] = 12
         
cluster_counts = cdata.obs.groupby(['condition', 'neuron_cluster']).size().reset_index(name='Count')
print(cluster_counts)

## Calculate the total number of neurons in each group
total_neurons_per_group = cdata.obs.groupby('condition').size().reset_index(name='Total_Neurons_group')
print(total_neurons_per_group)

## Merge the total neuron counts into the cluster counts DataFrame
cluster_counts = pd.merge(cluster_counts, total_neurons_per_group, on='condition')
cluster_counts

## Calculate the frequency of each neuron cluster within each group
cluster_counts['Frequency'] = cluster_counts['Count'] / cluster_counts['Total_Neurons_group']
cluster_counts
## Plot the cluster frequencies using seaborn
plt.figure(figsize=(3, 5))
g = sns.catplot(
    x='neuron_cluster', 
    y='Frequency',  # Use the calculated frequency
    hue='condition', 
    data=cluster_counts, 
    kind='bar', 
    height=3, 
    aspect=2, 
    palette='viridis'
)

# Customize the plot
g.set_axis_labels(" ", "% Neuronal Population\nper Group")  # Custom x and y labels
g._legend.set_title('Group')
g._legend.set_bbox_to_anchor((0.95, .85))  # Adjust the legend position
plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.savefig("./figures/cluster_frequencies.pdf", format='pdf', bbox_inches='tight') 
plt.show()



###############################################################################

## PATHOLOGY-ASSOCIATED GENES HEATMAP (FROM HOLTZMANN PAPER)

plt.rcParams['font.size'] = 16

## With Zbtb20
genes = ['Arpp21', 'R3hdm1', 'Rorb', 'Cux1', 'Cux2', 'Mef2c', 'Zbtb20']

### Six groups
sc.pl.matrixplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-six-groups-Zbtb20.pdf')

sc.pl.dotplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-six-groups-Zbtb20.pdf', cmap = 'viridis')

### Four groups
sc.pl.matrixplot(cdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-four-groups-Zbtb20.pdf')

sc.pl.matrixplot(cdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-four-groups-Zbtb20.pdf')


## Without Zbtb20
genes = ['Arpp21', 'R3hdm1', 'Rorb', 'Cux1', 'Cux2', 'Mef2c']

### Six groups
sc.pl.matrixplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-six-groups.pdf')

sc.pl.dotplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-six-groups.pdf', cmap = 'viridis')

### Four groups
sc.pl.matrixplot(cdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-four-groups.pdf')

sc.pl.matrixplot(cdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-four-groups.pdf')


###############################################################################

# NEURONAL HEALTH-ASSOCIATED GENES HEATMAP (HOLTZMANN PAPER)

genes = ['Prox1', 'Synpr', 'C1ql2', 'C1ql3', 'Camk2a', 'Camk2b', 'Tmem108', 'Ppfia2', 'Rfx3', 'Lrrtm4', 'Btbd9', 'Cntnap5a', 'Erc2']

plt.rcParams['font.size'] = 16


### Six groups
sc.pl.matrixplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'neuronal-health-genes-six-groups.pdf')

sc.pl.dotplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'neuronal-health-genes-six-groups.pdf', cmap = 'viridis')

### Four groups
sc.pl.matrixplot(cdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'neuronal-health-genes-four-groups.pdf')

sc.pl.dotplot(cdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'neuronal-health-genes-four-groups.pdf',  cmap = 'viridis')


###############################################################################

# GENES BY TREATMENT GROUP

genes = ['Prox1', 'Synpr', 'C1ql2', 'C1ql3', 'Camk2a', 'Camk2b', 'Tmem108', 'Ppfia2', 'Rfx3', 'Lrrtm4', 'Btbd9', 'Cntnap5a', 'Erc2']

plt.rcParams['font.size'] = 16


### Six groups
sc.pl.matrixplot(bdata, var_names = genes, groupby = 'condition', standard_scale = 'var', dendrogram = True,  save = 'neuronal-health-genes-six-groups-condition.pdf')

sc.pl.dotplot(bdata, var_names = genes, groupby = 'condition', standard_scale = 'var',  save = 'neuronal-health-genes-six-groups-condition.pdf', cmap = 'viridis')

### Four groups
sc.pl.matrixplot(cdata, var_names = genes, groupby = 'condition', standard_scale = 'var',  save = 'neuronal-health-genes-four-groups-condition.pdf')

sc.pl.dotplot(cdata, var_names = genes, groupby = 'condition', standard_scale = 'var', save = 'neuronal-health-genes-four-groups-condition.pdf',  cmap = 'viridis')


###############################################################################