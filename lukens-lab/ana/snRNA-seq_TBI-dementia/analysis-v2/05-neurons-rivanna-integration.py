#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 18:14:07 2024

@author: maureen
"""


import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import scvi 

plt.rcParams['figure.dpi'] = 200
plt.rcParams['font.family'] = 'Arial'

###############################################################################

# IMPORT

## AnnData
data_dir = "/scratch/mnc3ra/2024-06_tbi-snseq-rivanna/h5ad/"

adata = sc.read_h5ad(os.path.join(data_dir, '4-tbi-annotated-full.h5ad'))

adata.obs['cell_type'].replace({'Ttr+': 'Choroid-plexus epithelial'}, inplace = True)

###############################################################################

# PREPARE ANNDATA

## Create neuronal subset

groups = ['Excitatory neuron', 'Inhibitory neuron']
bdata = adata[adata.obs['cell_type'].isin(groups)].copy()
sc.pp.highly_variable_genes(bdata, n_top_genes = 3000, subset = True)


groups = ['Excitatory neuron', 'Inhibitory neuron']
cdata = adata[adata.obs['cell_type'].isin(groups)].copy()

groups = ['Excitatory neuron', 'Inhibitory neuron', 'Choroid-plexus epithelial']
ddata = adata[adata.obs['cell_type'].isin(groups)].copy()
sc.pp.highly_variable_genes(ddata, n_top_genes = 3000, subset = True)

###############################################################################

# SET UP AND TRAIN MODELS

## Model 1
scvi.model.SCVI.setup_anndata(
    bdata,
    layer='counts',
    categorical_covariate_keys=['group'],
    continuous_covariate_keys=['total_counts', 'pct_counts_mt', 'pct_counts_ribosomal'],
)

model_1 = scvi.model.SCVI(bdata)

scvi.train.Trainer(accelerator='cpu', devices=1)
model_1.train()


## Model 2
scvi.model.SCVI.setup_anndata(
    cdata,
    layer='counts',
    categorical_covariate_keys=['group'],
    continuous_covariate_keys=['total_counts', 'pct_counts_mt', 'pct_counts_ribosomal'],
)

model_2 = scvi.model.SCVI(cdata)

scvi.train.Trainer(accelerator='cpu', devices=1)
model_2.train()


## Model 3
scvi.model.SCVI.setup_anndata(
    ddata,
    layer='counts',
    categorical_covariate_keys=['group'],
    continuous_covariate_keys=['total_counts', 'pct_counts_mt', 'pct_counts_ribosomal'],
)

model_3 = scvi.model.SCVI(ddata)

scvi.train.Trainer(accelerator='cpu', devices=1)
model_3.train()



################################################################################################################################

# SAVE / IMPORT MODEL

scvi_dir = '/scratch/mnc3ra/2024-06_tbi-snseq-rivanna/scvi/'

## Save model_1
model_1_dir = os.path.join(scvi_dir, 'model_1')
print(model_1_dir)
#model_1.save(model_1_dir)


## Save model_2
model_2_dir = os.path.join(scvi_dir, 'model_2')
print(model_2_dir)
#model_2.save(model_2_dir)

## Save model_3
model_3_dir = os.path.join(scvi_dir, 'model_3')
print(model_3_dir)
#model_3.save(model_3_dir)

####


## Import model

model_1_dir = os.path.join(scvi_dir, 'model_1')
print(model_1_dir)
model_1 = scvi.model.SCVI.load(model_1_dir, adata=adata)
model_1

################################################################################################################################
