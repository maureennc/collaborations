#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 13:26:34 2024

@author: maureen
"""



import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

plt.rcParams['figure.dpi'] = 500
plt.rcParams['figure.figsize'] =(10, 5)

plt.rcParams['font.family'] = 'Arial'

###############################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '5-neuron-recluster-full.h5ad'))
sc.pl.umap(adata, color = 'neuron_cluster')

## Prepare data
adata = adata[adata.obs['neuron_cluster'] != 'Inhibitory 10'].copy()
adata.obs['neuron_cluster'].replace({'Mixed Ex/Inhib': 'Mixed Ex/In'}, inplace=True)

groups = ['A', 'B', 'C', 'E']
adata = adata[adata.obs['group'].isin(groups)].copy()
#adata.X = adata.layers['counts'].copy()
#sc.pp.normalize_total(adata)

###############################################################################

## HOLTZMAN PAPER EXCITATORY 10/"CLUSTER 6" GENES

genes = ['Arpp21', 'R3hdm1', 'Rorb', 'Cux1', 'Cux2', 'Mef2c'] #Zbtb20


sc.pl.matrixplot(adata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True)


sc.tl.dendrogram(adata, groupby = 'group')


excitatory = adata[adata.obs['cell_type'] == 'Excitatory neuron'].copy()
sc.pl.matrixplot(excitatory, var_names = genes, groupby = 'group', standard_scale = 'var', dendrogram = False, swap_axes = True)


ex10 = adata[adata.obs['neuron_cluster'] == 'Excitatory 10'].copy()
sc.pl.violin(ex10, keys = 'Rorb', groupby = 'group')

ex7 = adata[adata.obs['neuron_cluster'] == 'Excitatory 7'].copy()
sc.pl.violin(ex7, keys = 'Rorb', groupby = 'group')


ex12 = adata[adata.obs['neuron_cluster'] == 'Excitatory 12'].copy()
sc.pl.violin(ex12, keys = 'Rorb', groupby = 'group')

###############################################################################