#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 10:53:13 2024

@author: maureen

2024-09-20 Figures
"""

import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
#sc.settings.figdir = save_dir


plt.rcParams['figure.dpi'] = 500
#plt.rcParams['figure.figsize'] =(10, 5)
#plt.rcParams['font.size'] = 18

plt.rcParams['font.family'] = 'Arial'
plt.rcdefaults()

###############################################################################

# IMPORT AND PREPARE DATA

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '5-neuron-recluster-full.h5ad'))

## Prepare data
adata = adata[adata.obs['neuron_cluster'] != 'Inhibitory 10'].copy()
adata.obs['neuron_cluster'].replace({'Mixed Ex/Inhib': 'Mixed Ex/In'}, inplace=True)


###############################################################################

# HEATMAP 1 - HOLTZMANN GENES - ALL GROUPS

## 1) update the current one with the highly expressed genes in Ex 10 by adding Brinp3 and the repressor Zbtb20.

save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/figures/2024-09-20_figures'
sc.settings.figdir = save_dir


## Plotting
genes = ['Arpp21', 'R3hdm1', 'Rorb', 'Cux1', 'Cux2', 'Mef2c', 'Brinp3', 'Zbtb20'] 

plt.rcParams['font.size'] = 16
sc.pl.matrixplot(adata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, figsize = None)


plt.rcParams['font.size'] = 14

sc.pl.matrixplot(adata, var_names = genes, groupby = 'group', standard_scale = 'var', dendrogram = False, swap_axes = True, figsize = None)

###############################################################################

# HEATMAP 1 - HOLTZMANN GENES - PRIORITY GROUPS

groups = ['A', 'B', 'C', 'E']
bdata = adata[adata.obs['group'].isin(groups)].copy()

## Plotting
genes = ['Arpp21', 'R3hdm1', 'Rorb', 'Cux1', 'Cux2', 'Mef2c', 'Brinp3', 'Zbtb20'] 

plt.rcParams['font.size'] = 16
sc.pl.matrixplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, figsize = None, save = 'holtzmann-genes.pdf')

plt.rcParams['font.size'] = 16


bata = bdata[bdata.obs['neuron_cluster'] == 'Excitatory neuron 10'].copy()
## Prepare labels for group / condition
condition = {
    'A': 'Sham + GFP',
    'B': 'Sham + VEGFC',
    'C': 'TBI + GFP',
    'E': 'TBI + VEGFC',
}
bdata.obs['condition'] = bdata.obs['group'].map(condition)
sc.pl.matrixplot(bdata, var_names = genes, groupby = 'condition', standard_scale = 'var', dendrogram = False, swap_axes = True, save = 'holtzmann-genes-by-group_subset.pdf')


## ON EXCITATORY 10

sc.pl.matrixplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, figsize = None)



## Gene scoring

#sc.tl.score_genes(bdata, gene_list=genes, score_name='Pathology-associated Signature')
#data = pd.DataFrame(bdata.obs[['Pathology-associated Signature']])

#sc.pl.violin(bdata, keys='Pathology-associated Signature', groupby='neuron_cluster', jitter=0.4, rotation=90)

###############################################################################

# HEATMAP 2 - HEALTHY NEURONAL GENES

## Plotting
genes = ['Prox1', 'Synpr', 'C1ql2', 'C1ql3', 'Camk2a', 'Camk2b', 'Tmem108', 'Ppfia2', 'Rfx3', 'Lrrtm4', 'Btbd9', 'Cntnap5a', 'Erc2']

plt.rcParams['font.size'] = 16

sc.pl.matrixplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, figsize = None, save = 'neuronal-health-genes.pdf')


## Neuronal health genes by group

sc.pl.matrixplot(bdata, var_names = genes, groupby = 'condition', standard_scale = 'var', dendrogram = False, swap_axes = True, figsize = None, save = 'neuronal-health-genes-by-group.pdf')


sc.pl.dotplot(bdata, var_names = genes, groupby = 'condition', standard_scale = 'var', dendrogram = False, swap_axes = True)




###############################################################################

# GENES OF INTEREST

csv_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/cell_type/wald/pairwise-comps/size-factors'

de_results = pd.read_csv(os.path.join(csv_dir, 'wald-test-cell_type-size-factors-full.csv'))

comparisons = ['AC', 'CE']

genes = ['Arpp21', 'R3hdm1', 'Rorb', 'Cux1', 'Cux2', 'Mef2c', 'Brinp3', 'Zbtb20', 'Prox1', 'Synpr', 'C1ql2', 'C1ql3', 'Camk2a', 'Camk2b', 'Tmem108', 'Ppfia2', 'Rfx3', 'Lrrtm4', 'Btbd9', 'Cntnap5a', 'Erc2'] 

significant_hits = de_results[
    (de_results['comparison'].isin(comparisons)) &
    (de_results['gene'].isin(genes)) &
    (de_results['qval'] < 0.05) &
    (de_results['cell_type'] == 'Excitatory neuron')
].copy()

significant_hits.to_csv(os.path.join(csv_dir, 'gene-lists/excitatory-sig-genes-of-interest.csv'), index = False)

###############################################################################