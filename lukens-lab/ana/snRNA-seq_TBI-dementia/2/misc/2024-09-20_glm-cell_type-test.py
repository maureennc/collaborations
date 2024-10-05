#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 17:42:28 2024

@author: maureen
"""

import os
import scanpy as sc
from scipy.sparse import issparse
import diffxpy.api as de
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 200

###############################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '4-tbi-annotated-full.h5ad'))

###############################################################################

# ADD TESTING ANNOTATIONS

## TBI
tbi_mapping = {
    'A': 'Sham',
    'B': 'Sham',
    'C': 'TBI',
    'D': 'TBI',
    'E': 'TBI',
    'F': 'TBI'
}

## AAV
aav_mapping = {
    'A': 'GFP',
    'B': 'VEGFC',
    'C': 'GFP',
    'D': 'GFP',
    'E': 'VEGFC',
    'F': 'VEGFC'
}

## Ipsi-contra
side_mapping = {
    'A': 'Ipsilateral',
    'B': 'Ipsilateral',
    'C': 'Ipsilateral',
    'D': 'Contralateral',
    'E': 'Ipsilateral',
    'F': 'Contralateral'
}

adata.obs['TBI'] = adata.obs['group'].map(tbi_mapping)
adata.obs['AAV'] = adata.obs['group'].map(aav_mapping)
adata.obs['side'] = adata.obs['group'].map(side_mapping)

adata.obs[['group', 'TBI', 'AAV', 'side']].head()

###############################################################################

# PARAMETERIZED FUNCTIONS 

## Calculate size factors
def calculate_size_factors(subset_data):
    total_counts = subset_data.X.sum(axis=1)
    size_factors = total_counts / np.median(total_counts)
    return size_factors


skipped_analyses = []

## TBI
def run_de_tbi_effects_analysis(adata, cell_type, min_cells_per_group=30, max_imbalance_ratio=0.9):
    # Filter data for the specific cell type
    subset_data = adata[adata.obs['cell_type'] == cell_type].copy()

    # Checking for minimum cell count per group within the specific cell type
    group_counts = subset_data.obs['TBI'].value_counts()
    if any(group_counts < min_cells_per_group):
        skipped_analyses.append({
            'cell_type': cell_type, 'comparison': 'TBI', 'reason': 'Not enough cells in one or more groups'
        })
        return "Skipped: not enough cells in one or more groups for TBI analysis"

    # Checking for imbalance in cell distribution between groups
    if max(group_counts) / sum(group_counts) > max_imbalance_ratio:
        skipped_analyses.append({
            'cell_type': cell_type, 'comparison': 'TBI', 'reason': 'Imbalance between groups'
        })
        return "Skipped: imbalance between groups for TBI analysis"

    # Use raw counts
    subset_data.X = subset_data.layers['counts'].copy()

    # Filter genes expressed in at least 10 cells
    sc.pp.filter_genes(subset_data, min_cells=20)
    print(subset_data)

    if issparse(subset_data.X):
        subset_data.X = subset_data.X.toarray()

    ## Calculate size factors
    size_factors = calculate_size_factors(subset_data)
    subset_data.obs['size_factors'] = size_factors

    ## Run Wald test for TBI
    try:
        test_result = de.test.wald(
            data=subset_data,
            formula_loc="~ 1 + TBI + AAV + TBI:AAV",  # 1 + TBI + AAV led to convergence issues for neurons
            factor_loc_totest=["TBI"]
        )
    except Exception as e:
        return f"Error in DE test for TBI: {str(e)}"

    # Collect results
    df_result = test_result.summary()
    df_result['cell_type'] = cell_type
    df_result['comparison'] = "TBI"

    return df_result


## AAV
def run_de_aav_effects_analysis(adata, cell_type, min_cells_per_group=30, max_imbalance_ratio=0.9):
    # Filter data for the specific cell type
    subset_data = adata[adata.obs['cell_type'] == cell_type].copy()

    # Checking for minimum cell count per group within the specific cell type
    group_counts = subset_data.obs['AAV'].value_counts()
    if any(group_counts < min_cells_per_group):
        skipped_analyses.append({
            'cell_type': cell_type, 'comparison': 'AAV', 'reason': 'Not enough cells in one or more groups'
        })
        return "Skipped: not enough cells in one or more groups for AAV analysis"

    # Checking for imbalance in cell distribution between groups
    if max(group_counts) / sum(group_counts) > max_imbalance_ratio:
        skipped_analyses.append({
            'cell_type': cell_type, 'comparison': 'AAV', 'reason': 'Imbalance between groups'
        })
        return "Skipped: imbalance between groups for AAV analysis"

    # Use raw counts
    subset_data.X = subset_data.layers['counts'].copy()

    # Filter genes expressed in at least 10 cells
    sc.pp.filter_genes(subset_data, min_cells=20)
    print(subset_data)

    if issparse(subset_data.X):
        subset_data.X = subset_data.X.toarray()

    ## Calculate size factors
    size_factors = calculate_size_factors(subset_data)
    subset_data.obs['size_factors'] = size_factors

    ## Run Wald test for AAV
    try:
        test_result = de.test.wald(
            data=subset_data,
            formula_loc="~ 1 + AAV + TBI + TBI:AAV",
            factor_loc_totest=["AAV"]
        )
    except Exception as e:
        return f"Error in DE test for AAV: {str(e)}"

    # Collect results
    df_result = test_result.summary()
    df_result['neuron_cluster'] = cell_type
    df_result['comparison'] = "AAV"

    return df_result


## INTERACTION
def run_de_interaction_effects_analysis(adata, cell_type, min_cells_per_group=30, max_imbalance_ratio=0.9):
    # Filter data for the specific cell type
    subset_data = adata[adata.obs['cell_type'] == cell_type].copy()

    # Checking for minimum cell count per group within the specific cell type
    group_counts = subset_data.obs.groupby(['AAV', 'TBI']).size()
    if any(group_counts < min_cells_per_group):
        skipped_analyses.append({
            'cell_type': cell_type, 'comparison': 'TBI:AAV', 'reason': 'Not enough cells in one or more groups'
        })
        return "Skipped: not enough cells in one or more groups for interaction analysis"

    # Checking for imbalance in cell distribution between groups
    if max(group_counts) / sum(group_counts) > max_imbalance_ratio:
        skipped_analyses.append({
            'cell_type': cell_type, 'comparison': 'TBI:AAV', 'reason': 'Imbalance between groups'
        })
        return "Skipped: imbalance between groups for interaction analysis"

    # Use raw counts
    subset_data.X = subset_data.layers['counts'].copy()

    # Filter genes expressed in at least 10 cells
    sc.pp.filter_genes(subset_data, min_cells=20)
    print(subset_data)

    if issparse(subset_data.X):
        subset_data.X = subset_data.X.toarray()

    ## Calculate size factors
    size_factors = calculate_size_factors(subset_data)
    subset_data.obs['size_factors'] = size_factors

    ## Run Wald test for the interaction between AAV and TBI
    try:
        test_result = de.test.wald(
            data=subset_data,
            formula_loc="~ 1 + TBI + AAV + TBI:AAV",
            factor_loc_totest=["TBI:AAV"] 
        )
    except Exception as e:
        return f"Error in DE test for interaction: {str(e)}"

    # Collect results
    df_result = test_result.summary()
    df_result['cell_type'] = cell_type
    df_result['comparison'] = "TBI:AAV"

    return df_result

# After all analyses, convert skipped analyses to a DataFrame
def get_skipped_analyses():
    df_skipped = pd.DataFrame(skipped_analyses)
    return df_skipped


###############################################################################

# RUN DE FOR TBI, AAV, AND INTERACTION EFFECTS

cell_types = ['Excitatory neuron']

## TBI
all_results_tbi = []

for cell_type in cell_types:
    result_tbi = run_de_tbi_effects_analysis(adata=adata, cell_type=cell_type)
    if isinstance(result_tbi, pd.DataFrame):
        all_results_tbi.append(result_tbi)
    else:
        print(f"Skipping TBI analysis for {cell_type}: {result_tbi}")
        
results_tbi = pd.concat(all_results_tbi, ignore_index=True)


## AAV
all_results_aav = []

for cell_type in cell_types:
    result_aav = run_de_aav_effects_analysis(adata=adata, cell_type=cell_type)
    if isinstance(result_aav, pd.DataFrame):
        all_results_aav.append(result_aav)
    else:
        print(f"Skipping AAV analysis for {cell_type}: {result_aav}")

results_aav = pd.concat(all_results_aav, ignore_index=True)


## Interaction
all_results_interaction = []

for cell_type in cell_types:
    result_interaction = run_de_interaction_effects_analysis(adata=adata, cell_type=cell_type)
    if isinstance(result_interaction, pd.DataFrame):
        all_results_interaction.append(result_interaction)
    else:
        print(f"Skipping interaction analysis for {cell_type}: {result_interaction}")
        
results_interaction = pd.concat(all_results_interaction, ignore_index=True)


## Concatenate full results
final_results = pd.concat([results_tbi, results_aav, results_interaction], ignore_index=True)

## Convert skipped_analysis to dataframe
skip = pd.DataFrame(skipped_analyses)

###############################################################################

# FILTER RESULTS

def filter_de_results(df, max_log2fc=100):
    return df[(abs(df['log2fc']) <= max_log2fc)]

filtered_results_tbi = filter_de_results(results_tbi)
filtered_results_aav = filter_de_results(results_aav)
filtered_results_interaction = filter_de_results(results_interaction)

filtered_tbi = filtered_results_tbi[filtered_results_tbi['mean'] > 0.1].copy()
###############################################################################

# EXPORT DE RESULTS

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/multivariate/excitatory/1"

results_aav.to_csv(os.path.join(save_dir, 'aav-main-effect-excitatory-full.csv'), index = False)
filtered_results_aav.to_csv(os.path.join(save_dir, 'aav-main-effect-excitatory-filtered.csv'), index = False)

results_tbi.to_csv(os.path.join(save_dir, 'tbi-main-effect-excitatory-full.csv'), index = False)
filtered_results_tbi.to_csv(os.path.join(save_dir, 'tbi-main-effect-excitatory-filtered.csv'), index = False)

results_interaction.to_csv(os.path.join(save_dir, 'interaction-effect-excitatory-full.csv'), index = False)
filtered_results_interaction.to_csv(os.path.join(save_dir, 'interaction-effect-excitatory-filtered.csv'), index = False)

#skip.to_csv(os.path.join(save_dir, 'log.csv'), index = False)

###############################################################################
