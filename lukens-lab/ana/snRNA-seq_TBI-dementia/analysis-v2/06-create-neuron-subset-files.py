#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 10:07:53 2024

@author: maureen
"""

import os
import pandas as pd
import scanpy as sc
import scvi
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['figure.dpi'] = 200
plt.rcParams['font.family'] = 'Arial'

###############################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '4-tbi-annotated-full.h5ad'))
bdata = sc.read_h5ad(os.path.join(data_dir, '4-tbi-annotated-full.h5ad'))

###############################################################################

# PREPARE ANNDATA

## Create subset
groups = ['Excitatory neuron', 'Inhibitory neuron']
adata = adata[adata.obs['cell_type'].isin(groups)].copy()
sc.pp.highly_variable_genes(adata, n_top_genes = 3000, subset = True)

###############################################################################

# PROCESS ANNDATA

scvi_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/scvi"

## Import scVI model
model_dir = os.path.join(scvi_dir, 'model_1')
print(model_dir)
model = scvi.model.SCVI.load(model_dir, adata=adata)
model

## Add scvi latent key to obsm
SCVI_LATENT_KEY = "X_scVI"

## For model 1 (neuronal HVGs)
latent = model.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
print(f"Model latent shape: {latent.shape}")


###############################################################################

# INITIAL CLUSTERING

sc.pp.neighbors(adata, use_rep = 'X_scVI', random_state = 0)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added = 'leiden_scVI2', resolution = 1)
sc.pl.umap(adata, color = 'leiden_scVI2', legend_loc = 'on data')

genes = ['Rbfox3', 'Map2', 'Syn1', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 
         'Cx3cr1', 'Hexb', 'Mrc1',
         'Slc1a2', 'Slc1a3', 'Gfap', 'Aldh1l1', 'Aqp4', 'Atp1a2',
         'Pecam1', 'Tie1', 'Cldn5', 'Flt1',
         'Pdgfrb', 'Rgs5', 
         'Pdgfra', 'Cspg4', 'Sox10',
         'Mbp', 'Plp1', 'Mog'] 

sc.pl.dotplot(adata, genes, groupby = 'leiden_scVI2', standard_scale = 'var', use_raw = True, dendrogram = True)


###############################################################################

# ANNOTATION

neuron_cluster = {str(i): "" for i in range(30)}
neuron_cluster

neuron_cluster = {
 '0': 'Excitatory 1',
 '1': 'Excitatory 2',
 '2': 'Excitatory 3',
 '3': 'Inhibitory 1',
 '4': 'Excitatory 4',
 '5': 'Mixed Ex/Inhib', # Formerly labeled "Inhibitory 2"
 '6': 'Inhibitory 2',
 '7': 'Inhibitory 3',
 '8': 'Excitatory 5',
 '9': 'Excitatory 6',
 '10': 'Inhibitory 4',
 '11': 'Excitatory 7',
 '12': 'Excitatory 8',
 '13': 'Excitatory 9',
 '14': 'Inhibitory 5',
 '15': 'Excitatory 10',
 '16': 'Inhibitory 6',
 '17': 'Inhibitory 7',
 '18': 'Inhibitory 8',
 '19': 'Excitatory 11',
 '20': 'Excitatory 12',
 '21': 'Excitatory 13',
 '22': 'Excitatory 14',
 '23': 'Inhibitory 9',
 '24': 'Excitatory 15',
 '25': 'Excitatory 16',
 '26': 'Excitatory 17',
 '27': 'Excitatory 18',
 '28': 'Excitatory 19',
 '29': 'Inhibitory 10'}

adata.obs['neuron_cluster'] = adata.obs.leiden_scVI2.map(neuron_cluster)

sc.pl.umap(adata, color = 'neuron_cluster')


###############################################################################

# MAPPING - FIND INTERSECTION OF BARCODES

common_barcodes = adata.obs_names.intersection(bdata.obs_names)
bdata = bdata[common_barcodes].copy()

###############################################################################

# MAPPING - TRANSFER METADATA TO GENOME-SCALE BRANCH

## adata.obs
unique_columns = [col for col in adata.obs.columns if col not in bdata.obs.columns]

bdata.obs = bdata.obs.join(adata.obs[unique_columns], how='left')

## uns entries
for key in adata.uns.keys():
    #if key not in bdata.uns:
        bdata.uns[key] = adata.uns[key]

## obsm entries
for key in adata.obsm.keys():
    #if key not in bdata.obsm:
        bdata.obsm[key] = adata.obsm[key]

## obsp entries
for key in adata.obsp.keys():
    #if key not in bdata.obsp:
        bdata.obsp[key] = adata.obsp[key]

###############################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/h5ad"

adata.X = csr_matrix(adata.X)
bdata.X = csr_matrix(bdata.X)

adata.write_h5ad(os.path.join(save_dir, '5-neuron-recluster-hvg.h5ad'))
bdata.write_h5ad(os.path.join(save_dir, '5-neuron-recluster-full.h5ad'))

###############################################################################