#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:00:22 2024

@author: maureen
"""


import os
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import diffxpy.api as de
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

plt.rcParams['figure.dpi'] = 200
plt.rcParams['font.family'] = 'Arial'

###############################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '4-tbi-annotated-full.h5ad'))

###############################################################################

# PREPARE ANNDATA

## Mask artifact clusters
adata = adata[adata.obs['cell_type'] != 'Unknown'].copy()

## Rename Ttr+
adata.obs['cell_type'].replace({'Ttr+': 'Choroid-plexus epithelial'}, inplace = True)

## Restore the raw counts for updated normalization
adata.X = adata.layers['counts'].copy()

## Convert to dense array
type(adata.X)
adata.X = adata.X.toarray()
type(adata.X)

###############################################################################

# PARAMETERIZED FUNCTION

## size factors function
def calculate_size_factors(subset_data):
    total_counts = subset_data.X.sum(axis=1)
    size_factors = total_counts / np.mean(total_counts)
    return size_factors

## DE function
def run_de_analysis(adata, cell_type, reference_group, comparison_group, min_cells_per_group=40, max_imbalance_ratio=0.9):
    # Filter data for the specific cell type and groups
    subset_data = adata[(adata.obs['cell_type'] == cell_type) & 
                        (adata.obs['group'].isin([reference_group, comparison_group]))].copy()

    # Checking for minimum cell count per group within the specified cell type
    group_counts = subset_data.obs['group'].value_counts()
    if any(group_counts < min_cells_per_group):
        return f"Skipped: {reference_group} vs {comparison_group} (n is too low)"

    # Checking for imbalance in cell distribution between groups within the subset
    if max(group_counts) / sum(group_counts) > max_imbalance_ratio:
        return f"Skipped: {reference_group} vs {comparison_group} (imbalance)"

    # Proceed if checks pass
    subset_data.X = subset_data.layers['counts'].copy()

    # Filter genes that are expressed in at least 3 cells
    sc.pp.filter_genes(subset_data, min_cells=3)

    if issparse(subset_data.X):
        subset_data.X = subset_data.X.toarray()

    # Calculate size factors for each cell
    size_factors = calculate_size_factors(subset_data)
    
    # Add size factors to the 'obs' dataframe in adata
    subset_data.obs['size_factors'] = size_factors

    # Run the Wald test, passing the size factors directly
    test_result = de.test.wald(
        data=subset_data,
        formula_loc="~ 1 + group",
        factor_loc_totest="group",
        size_factors=size_factors  # Pass the size factors here
    )
    
    # Collect results
    df_result = test_result.summary()
    df_result['cell_type'] = cell_type
    df_result['comparison'] = f"{reference_group}{comparison_group}"
    
    return df_result


###############################################################################

# RUN DE

cell_types = ['Excitatory neuron', 'Inhibitory neuron', 'Astrocyte', 'Oligodendrocyte', 'OPC', 'Microglia', 'Choroid-plexus epithelial', 'Fibroblast']

#comparisons = [('A', 'C'), ('C', 'E'), ('A', 'B')]

comparisons = [('A', 'C'), ('C', 'E'), ('A', 'B'), ('A', 'D'), ('C', 'D'), ('D', 'F'), ('E', 'F')]

all_results = []

for cell_type in cell_types:
    for reference_group, comparison_group in comparisons:
        result = run_de_analysis(
            adata=adata,
            cell_type=cell_type,
            reference_group=reference_group,
            comparison_group=comparison_group
        )
        if isinstance(result, pd.DataFrame):
            all_results.append(result)
        else:
            print(f"Skipping {cell_type} {reference_group} vs {comparison_group}: {result}")  

# Concatenate results
final_results = pd.concat(all_results, ignore_index=True)

###############################################################################

# FILTER RESULTS

def filter_de_results(df, max_log2fc=100):
    return df[(abs(df['log2fc']) <= max_log2fc)]


filtered_results = filter_de_results(final_results)

###############################################################################

# EXPORT DE RESULTS

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/cell_type/wald/pairwise-comps/size-factors"

final_results.to_csv(os.path.join(save_dir, 'wald-test-cell_type-size-factors-full.csv'), index = False)

filtered_results.to_csv(os.path.join(save_dir, 'wald-test-cell_type-size-factors-filtered.csv'), index = False)

###############################################################################

###############################################################################

# IMPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/cell_type/wald/pairwise-comps/size-factors"


final_results =  pd.read_csv(os.path.join(save_dir, 'wald-test-cell_type-size-factors-full.csv'))

filtered_results =  pd.read_csv(os.path.join(save_dir, 'wald-test-cell_type-size-factors-filtered.csv'))


save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/wald/pairwise-comps/size-factors"

## Export filtered versions for volcano plots
comparison_ids = filtered_results['comparison'].unique()

# Loop through each comparison and save the filtered data
for comp in comparison_ids:
    # Filter the DataFrame for the current comparison
    df_filtered = filtered_results[filtered_results['comparison'] == comp]
    
    # Define the file name using the comparison ID
    file_name = f"{comp}-comparison-cell_type-filtered.csv"
    
    # Save the DataFrame as a CSV file
    #df_filtered.to_csv(os.path.join(save_dir, file_name), index=False)

    print(f"Saved: {file_name}")

###############################################################################

# VISUALIZE

from functools import reduce

# List of comparisons
comparisons = ['AC', 'CE', 'AB']
#comparisons = ['AC', 'CE', 'AB', 'AD', 'CD', 'DF', 'EF']

# Specify order of cell types for heatmap
order = ['Astrocyte', 'Excitatory neuron', 'Inhibitory neuron', 'Microglia', 'OPC', 'Oligodendrocyte', 'Fibroblast', 'Choroid-plexus epithelial' ]

# Create a reference DataFrame for cell types
cell_type_reference = pd.DataFrame(order, columns=['cell_type'])

# Initialize an empty list to store DEG counts per comparison
data_list = []

# Define the mean expression threshold
mean_expression_threshold = 0.1

# Loop over each comparison
for comp in comparisons:
    # Filter for the specific comparison in filtered_results
    df_comp = filtered_results[(filtered_results['comparison'] == comp)]

    # Additional filter to apply the mean expression threshold
    df_comp = df_comp[df_comp['mean'] > mean_expression_threshold]

    # Filter for significant DEGs with p-value < 0.05 and abs(log2FC) > 0.5
    df_significant = df_comp[(df_comp['qval'] < 0.05) & (abs(df_comp['log2fc']) > 0.5)]
    
    # Count the number of significant DEGs per cell type
    deg_counts = df_significant.groupby('cell_type').size().reset_index(name=f'{comp}')
    
    # Merge with the cell type reference to ensure all cell types are present
    deg_counts = pd.merge(cell_type_reference, deg_counts, on='cell_type', how='left').fillna(0)
    
    # Append to the list
    data_list.append(deg_counts)


# Merge all comparison dataframes on 'cell_type'
combined_df = reduce(lambda left, right: pd.merge(left, right, on='cell_type', how='outer'), data_list)
combined_df.fillna(0, inplace=True)  # Replace NaN with 0 where there are no DEGs

# Reorder the dataframe according to the desired cell type order
combined_df['cell_type'] = pd.Categorical(combined_df['cell_type'], categories=order, ordered=True)
combined_df = combined_df.sort_values('cell_type')

# Plot the heatmap
plt.figure(figsize=(9,6))
sns.heatmap(combined_df.set_index('cell_type'), annot=True, cmap='viridis', fmt="g")
plt.title('')
plt.ylabel(' ')
plt.xlabel('')
plt.grid(False)
plt.tight_layout()
plt.show()


###############################################################################

microglia = filtered_results[filtered_results['cell_type'] == "Microglia"].copy()
microglia = microglia[microglia['comparison'] == "AC"].copy()
microglia = microglia[microglia['qval'] < 0.05].copy()
microglia = microglia[microglia['mean'] > 0.25].copy()


inhibitory = filtered_results[filtered_results['cell_type'] == "Inhibitory neuron"].copy()
inhibitory = inhibitory[inhibitory['comparison'] == "CE"].copy()
inhibitory = inhibitory[inhibitory['qval'] < 0.05].copy()


excitatory = filtered_results[filtered_results['cell_type'] == "Excitatory neuron"].copy()
excitatory = excitatory[excitatory['comparison'] == "CE"].copy()
excitatory = excitatory[excitatory['qval'] < 0.05].copy()
