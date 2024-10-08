#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 11:36:03 2024

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


###############################################################################

# SETTINGS

os.chdir('/Users/maureen/Documents/projects/lukens-lab/nick/syk-aging-cosmx/analysis/4')

###############################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/nick/syk-aging-cosmx/analysis/3/h5ad/"

adata = sc.read_h5ad(os.path.join(data_dir, '2-syk-integrated-annotated.h5ad'))

###############################################################################

# CREATE DATA SUBSETS

bdata = adata[adata.obs['region'] == 'Hippocampus'].copy()
olig_lin = bdata[bdata.obs['cell_type'].isin(['Oligodendrocyte', 'OPC'])].copy()
oligo = bdata[bdata.obs['cell_type'] == 'Oligodendrocyte'].copy()
opc = bdata[bdata.obs['cell_type'] == 'OPC'].copy()
mg = bdata[bdata.obs['cell_type'] == 'Microglia'].copy()

###############################################################################

# AGE AND DISEASE-ASSOCIATED MICROGLIA AND OLIGO DOTPLOTS

oligo_genes = ['C4b', 'Hif3a', 'Cdkn1a', 'Sult1a1', 'Il12rb1', 'Stat1', 'Irf7', 'Irf9'] # removed Serpina3n

## Oligo-lineage
sc.pl.dotplot(olig_lin, oligo_genes, groupby = 'condition',  use_raw = True, standard_scale = 'var', save = 'oligo-lineage-standardized-expression-v2.pdf')
sc.pl.dotplot(olig_lin, oligo_genes, groupby = 'condition',  use_raw = True, dendrogram = False, save = 'oligo-lineage-log1p-expression-v2.pdf')
#sc.pl.dotplot(olig_lin, oligo_genes, groupby = 'condition',  use_raw = True, dendrogram = True, save = 'oligo-lineage-log1p-dendrogram.pdf')


## Mature Oligodendrocytes
sc.pl.dotplot(oligo, oligo_genes, groupby = 'condition',  use_raw = True, standard_scale = 'var', save = 'mature-oligo-standardized-expression-v2.pdf')
sc.pl.dotplot(oligo, oligo_genes, groupby = 'condition',  use_raw = True, dendrogram = False, save = 'mature-oligo-log1p-expression-v2.pdf')
#sc.pl.dotplot(oligo, oligo_genes, groupby = 'condition',  use_raw = True, dendrogram = 'True', save = 'mature-oligo-log1p-dendrogram.pdf')

## OPCs
sc.pl.dotplot(opc, oligo_genes, groupby = 'condition',  use_raw = True, standard_scale = 'var', save = 'opc-standardized-expression-v2.pdf')
sc.pl.dotplot(opc, oligo_genes, groupby = 'condition',  use_raw = True, dendrogram = False, save = 'opc-log1p-expression-v2.pdf')
#sc.pl.dotplot(opc, oligo_genes, groupby = 'condition',  use_raw = True, dendrogram = True, save = 'opc-log1p-dendrogram.pdf')


mg_genes = ['P2ry12', 'Cx3cr1', 'Tmem119', 'Cd33', 'Trem2', 'Tyrobp', 'Itgax', 'Cst7', 'Axl', 'Clec7a', 'Ctsb', 'Ctsd', 'Lpl', 'Lgals3bp', 'Syk', 'Lag3', 'Egr1'] # Moved Cx3cr1 up to homeostatic genes

## Microglia
sc.pl.dotplot(mg, mg_genes, groupby = 'condition',  use_raw = True, standard_scale = 'var', save = 'microglia-standardized-expression-v2.pdf')
sc.pl.dotplot(mg, mg_genes, groupby = 'condition',  use_raw = True, dendrogram = False, save = 'microglia-log1p-expression-v2.pdf')
#sc.pl.dotplot(mg, mg_genes, groupby = 'condition',  use_raw = True, dendrogram = True, save = 'microglia-log1p-dendrogram.pdf')

###############################################################################
