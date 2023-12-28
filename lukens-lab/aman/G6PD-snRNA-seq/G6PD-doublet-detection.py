#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 22:26:52 2023

@author: maureen
"""


import os
os.chdir("/Users/maureen/Documents/projects/lukens-lab/Aman/snRNAseq_G6PD_Project")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
sc.set_figure_params(scanpy = True, dpi = 150, dpi_save = 400)

import scrublet as scr
import anndata



##############################################

# IMPORT DATA

A1 = sc.read_10x_h5("h5_all/A1_filtered_feature_bc_matrix.h5")
A2 = sc.read_10x_h5("h5_all/A2_filtered_feature_bc_matrix.h5")
A3 = sc.read_10x_h5("h5_all/A3_filtered_feature_bc_matrix.h5")


A4 = sc.read_10x_h5("h5_5xFAD_Only/A4_filtered_feature_bc_matrix.h5")
A5 = sc.read_10x_h5("h5_5xFAD_Only/A5_filtered_feature_bc_matrix.h5")
A6 = sc.read_10x_h5("h5_5xFAD_Only/A6_filtered_feature_bc_matrix.h5")


A1.var_names_make_unique()
A2.var_names_make_unique()
A3.var_names_make_unique()

A4.var_names_make_unique()
A5.var_names_make_unique()
A6.var_names_make_unique()


#############################################


adata = A3
#change to adata file

##############################################


# DOUBLET REMOVAL

scrub = scr.Scrublet(adata.X, expected_doublet_rate = 0.20)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

scrub.plot_histogram()

#2D embedding
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))                                                                                  
scrub.plot_embedding('UMAP', order_points=True);
                
adata.obs['doublet_scores'] = doublet_scores
adata.obs['predicted_doublets'] = predicted_doublets

num_doublets = (adata.obs['doublet_scores'] > 0.5).sum()
print("Number of cells with doublet score > 0.5:", num_doublets)


##############################################

# WRITE

adata.write('A3_scrublet_anno.h5ad')


#adata.to_df().to_csv('expression_matrix.csv')
adata.obs.to_csv('A3_cell_metadata.csv')
#adata.var.to_csv('gene_metadata.csv')




###################################################


#expression_matrix = adata.to_df().T
# Save the expression matrix to CSV
#expression_matrix.to_csv('expression_matrix.csv')

# Cell Metadata
# Ensure the first column contains cell identifiers
#cell_metadata = adata.obs
#cell_metadata.insert(0, 'CellID', cell_metadata.index)
#cell_metadata.to_csv('cell_metadata.csv', index=False)

# Gene Metadata
# Ensure the first column contains gene identifiers
#gene_metadata = adata.var
#gene_metadata.insert(0, 'GeneID', gene_metadata.index)
#gene_metadata.to_csv('gene_metadata.csv', index=False)
####

#x = pd.read_csv('./deliverables/')


#y = x['predicted_doublets'] == True

#y.value_counts()



