#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:48:42 2024

@author: maureen
"""

import os
import pandas as pd
import scanpy as sc
import scvi
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['figure.dpi'] = 600
plt.rcParams['figure.figsize'] =(4, 3)

plt.rcParams['font.family'] = 'Arial'

###############################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '5-neuron-recluster-hvg.h5ad'))
adata

###############################################################################

# PROCESS ANNDATA

scvi_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/scvi"

## Import scVI model
model_dir = os.path.join(scvi_dir, 'model_1') # HVGs
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
sc.pl.umap(adata, color = 'group')


genes = ['Rbfox3', 'Map2', 'Syn1', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2'] 

sc.pl.dotplot(adata, genes, groupby = 'neuron_cluster', standard_scale = 'var', use_raw = True, dendrogram = True)


###############################################################################

# FUNCTION FOR DE RESULTS - PAIRWISE
save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/de-scvi/neuron'

comparisons = ['AC', 'CE', 'AB']

# Define the cell types (unique neuron clusters)
cell_types = adata.obs['neuron_cluster'].unique().tolist()

def run_scvi_de(cell_type, group1, group2, model):
    # Boolean indices for the two groups
    idx1 = (adata.obs['neuron_cluster'] == cell_type) & (adata.obs['group'] == group1)
    idx2 = (adata.obs['neuron_cluster'] == cell_type) & (adata.obs['group'] == group2)
    
    # Check if either group has fewer than 40 cells
    if idx1.sum() < 40 or idx2.sum() < 40:
        print(f"Skipping DE for {cell_type} between {group1} and {group2}: Not enough cells (<40).")
        return None
    
    # Check if either index is empty (i.e., no cells in that group)
    if idx1.sum() == 0 or idx2.sum() == 0:
        print(f"Skipping DE for {cell_type} between {group1} and {group2}: One of the groups is empty.")
        return None
    
    # Run differential expression if both groups have sufficient cells
    result = model.differential_expression(idx1=idx1, idx2=idx2)
    result['Cell_Type'] = cell_type  # Add a column for the cell type
    result.reset_index(inplace=True)  # Reset index to bring gene names into a column
    return result

## Collect results
for group1, group2 in comparisons:
    comparison_results = []
    for cell_type in cell_types:
        de_results = run_scvi_de(cell_type, group1, group2, model)
        if de_results is not None:  # Only append results if they exist
            comparison_results.append(de_results)
    
    # Concatenate results for all cell types for the current comparison
    if comparison_results:
        combined_results = pd.concat(comparison_results, ignore_index=True)
        file_name = f'{group1}{group2}-comparison.csv'
        full_path = os.path.join(save_dir, file_name)
        combined_results.to_csv(full_path, index=False)  # Ensure index=False is set if the gene names are already a column
        print(f"Saved: {full_path}")
    else:
        print(f"No valid DE results for {group1} vs {group2}")

###############################################################################

# PLOTTING

from functools import reduce
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Set result directory
results_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/de-scvi/neuron'

# Comparisons to be analyzed
comparisons = ['AC', 'CE', 'AB']

# Collect all unique cell types from all files
all_cell_types = set()
for comp in comparisons:
    file_path = os.path.join(results_dir, f"{comp}-comparison.csv")
    df = pd.read_csv(file_path)
    all_cell_types.update(df['Cell_Type'].unique())

all_cell_types = list(all_cell_types)  # Convert to list

# Prepare data for heatmap
data_list = []

for comp in comparisons:
    file_path = os.path.join(results_dir, f"{comp}-comparison.csv")
    df = pd.read_csv(file_path)
    print(f"Loaded {comp} with shape {df.shape}")

    # Filter significant DEGs based on FDR threshold, log2 fold change, and mean expression
    df_significant = df[
        (df['is_de_fdr_0.05'] == True) &  # Statistically significant
        (abs(df['lfc_median']) > 0.5) &  # Biologically relevant fold change
        ((df['raw_normalized_mean1'] > 0.01) | (df['raw_normalized_mean2'] > 0.01))  # Minimum mean expression threshold
    ]
    
    print(f"Filtered {comp} significant DEGs with shape {df_significant.shape}")

    # Group by cell type and count significant DEGs
    deg_counts = df_significant.groupby('Cell_Type').size().reset_index(name=f'{comp}')
    
    # Initialize a DataFrame for all cell types with zeros
    full_deg_counts = pd.DataFrame({'Cell_Type': all_cell_types})
    full_deg_counts = full_deg_counts.merge(deg_counts, on='Cell_Type', how='left').fillna(0)
    
    data_list.append(full_deg_counts)

# Merge all DataFrames on 'Cell_Type'
final_df = reduce(lambda left, right: pd.merge(left, right, on='Cell_Type', how='outer'), data_list)
final_df.fillna(0, inplace=True)

# Plot the heatmap
plt.figure(figsize=(5, 10))
sns.heatmap(final_df.set_index('Cell_Type'), annot=True, cmap='viridis', fmt='g')
plt.title('')
plt.ylabel('Cluster')
plt.xlabel('Comparison')
plt.grid(False)
plt.tight_layout()
plt.show()



###############################################################################

# WHITEOUT INSTEAD OF 0

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from functools import reduce
import os

results_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/de-scvi/neuron'
comparisons = ['AC', 'CE', 'AB']  # Example comparisons

# Collect all unique cell types from all files
all_cell_types = set()
for comp in comparisons:
    file_path = os.path.join(results_dir, f"{comp}-comparison.csv")
    df = pd.read_csv(file_path)
    all_cell_types.update(df['Cell_Type'].unique())

all_cell_types = list(all_cell_types)  # Convert to list

# Prepare data for heatmap
data_list = []

for comp in comparisons:
    file_path = os.path.join(results_dir, f"{comp}-comparison.csv")
    df = pd.read_csv(file_path)
    print(f"Loaded {comp} with shape {df.shape}")

    # Filter significant DEGs based on FDR threshold and log2 fold change
    df_significant = df[(df['is_de_fdr_0.05'] == True) & (abs(df['lfc_median']) > 0.5)]
    #df_significant = df_significant[(df['raw_normalized_mean1'] > 0.01) | (df['raw_normalized_mean2'] > 0.01)]

    print(f"Filtered {comp} significant DEGs with shape {df_significant.shape}")

    # Group by cell type and count significant DEGs
    deg_counts = df_significant.groupby('Cell_Type').size().reset_index(name=f'{comp}')
    
    # Initialize a DataFrame for all cell types, keeping NaN for missing values
    full_deg_counts = pd.DataFrame({'Cell_Type': all_cell_types})
    full_deg_counts = full_deg_counts.merge(deg_counts, on='Cell_Type', how='left')
    
    data_list.append(full_deg_counts)

# Merge all DataFrames on 'Cell_Type'
final_df = reduce(lambda left, right: pd.merge(left, right, on='Cell_Type', how='outer'), data_list)

# Plot the heatmap without filling NaNs (NaNs will appear as blank cells in the heatmap)
plt.figure(figsize=(9, 10))
sns.heatmap(final_df.set_index('Cell_Type'), annot=True, cmap='viridis', fmt='g', cbar=True)
plt.title('Significant DEGs per Cell Type Across Comparisons')
plt.ylabel('Cluster')
plt.xlabel('Comparison')
plt.grid(False)
plt.tight_layout()
plt.show()

###############################################################################

# DROP NaAN ROWS

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from functools import reduce
import os

# Set up paths and comparisons
results_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/de-scvi/neuron'
comparisons = ['AC', 'CE', 'AB']  # Example comparisons

# Collect all unique cell types from all files
all_cell_types = set()
for comp in comparisons:
    file_path = os.path.join(results_dir, f"{comp}-comparison.csv")
    df = pd.read_csv(file_path)
    all_cell_types.update(df['Cell_Type'].unique())

all_cell_types = list(all_cell_types)  # Convert to list

# Prepare data for heatmap
data_list = []

for comp in comparisons:
    file_path = os.path.join(results_dir, f"{comp}-comparison.csv")
    df = pd.read_csv(file_path)
    print(f"Loaded {comp} with shape {df.shape}")

    # Filter significant DEGs based on FDR threshold and log2 fold change
    df_significant = df[(df['is_de_fdr_0.05'] == True) & (abs(df['lfc_median']) > 0.5)]
    #df_significant = df_significant[(df['raw_normalized_mean1'] > 0.01) | (df['raw_normalized_mean2'] > 0.01)]

    print(f"Filtered {comp} significant DEGs with shape {df_significant.shape}")

    # Group by cell type and count significant DEGs
    deg_counts = df_significant.groupby('Cell_Type').size().reset_index(name=f'{comp}')
    
    # Initialize a DataFrame for all cell types, keeping NaN for missing values
    full_deg_counts = pd.DataFrame({'Cell_Type': all_cell_types})
    full_deg_counts = full_deg_counts.merge(deg_counts, on='Cell_Type', how='left')
    
    data_list.append(full_deg_counts)

# Merge all DataFrames on 'Cell_Type'
final_df = reduce(lambda left, right: pd.merge(left, right, on='Cell_Type', how='outer'), data_list)

# Drop rows with any NaN values
final_df.dropna(inplace=True)

# Plot the heatmap without NaNs (as these rows were dropped)
plt.figure(figsize=(5, 5))
sns.heatmap(final_df.set_index('Cell_Type'), annot=True, cmap='viridis', fmt='g', cbar=True)
plt.title('Significant DEGs per Cell Type Across Comparisons')
plt.ylabel('Cluster')
plt.xlabel('Comparison')
plt.grid(False)
plt.tight_layout()
plt.show()
