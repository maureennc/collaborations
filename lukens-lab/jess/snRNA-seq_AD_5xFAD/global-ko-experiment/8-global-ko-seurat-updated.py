#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 26 15:31:35 2024

@author: maureen
"""

import os
from pathlib import Path
import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse

################################################################################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/global-sting-ko/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, 'global-ko-trained-annotated.h5ad'))

################################################################################################################################

# PREPARE DATA

## Update group annotation
adata.obs['group'] = adata.obs['group'].replace('Sting-KO', 'Sting1-KO')

## Convert sparse matrix to dense
adata.X = adata.X.toarray()

################################################################################################################################

# EXPORT SEURATOBJECT COMPONENTS

save_dir = '/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/global-sting-ko/export/csv'

## Extract counts matrix
if isinstance(adata.X, scipy.sparse.spmatrix):
    counts = pd.DataFrame.sparse.from_spmatrix(adata.X, index=adata.obs_names, columns=adata.var_names)
else:
    counts = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)


## Export
counts.to_csv(os.path.join(save_dir, 'counts-matrix.csv'), sep=',')

adata.obs.to_csv(os.path.join(save_dir, 'metadata.csv'))
adata.var.to_csv(os.path.join(save_dir, 'features.csv'))

latent_representation = pd.DataFrame(adata.obsm['X_scVI'], index=adata.obs_names)
latent_representation.to_csv(os.path.join(save_dir, 'latent_representation.csv'))

################################################################################################################################

# EXPORT UPDATED H5AD

save_dir = '/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/global-sting-ko/export/h5ad'

adata.write_h5ad(os.path.join(save_dir, 'global-ko-processed-data.h5ad'))


################################################################################################################################

