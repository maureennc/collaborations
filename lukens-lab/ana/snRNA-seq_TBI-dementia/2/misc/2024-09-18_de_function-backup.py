#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 11:44:00 2024

@author: maureen
"""



## TBI
def run_de_tbi_effects_analysis(adata, cell_type, min_cells_per_group=40, max_imbalance_ratio=0.9):
    # Filter data for the specific cell type
    subset_data = adata[adata.obs['neuron_cluster'] == cell_type].copy()

    # Checking for minimum cell count per group within the specific cell type
    group_counts = subset_data.obs['TBI'].value_counts()
    if any(group_counts < min_cells_per_group):
        return f"Skipped: not enough cells in one or more groups for TBI analysis"

    # Checking for imbalance in cell distribution between groups
    if max(group_counts) / sum(group_counts) > max_imbalance_ratio:
        return f"Skipped: imbalance between groups for TBI analysis"

    # Use raw counts
    subset_data.X = subset_data.layers['counts'].copy()

    # Filter genes expressed in at least 10 cells
    sc.pp.filter_genes(subset_data, min_cells=10)
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
    df_result['neuron_cluster'] = cell_type
    df_result['comparison'] = "TBI"

    return df_result


## AAV
def run_de_aav_effects_analysis(adata, cell_type, min_cells_per_group=40, max_imbalance_ratio=0.9):
    # Filter data for the specific cell type
    subset_data = adata[adata.obs['neuron_cluster'] == cell_type].copy()

    # Checking for minimum cell count per group within the specific cell type
    group_counts = subset_data.obs['AAV'].value_counts()
    if any(group_counts < min_cells_per_group):
        return f"Skipped: not enough cells in one or more groups for AAV analysis"

    # Checking for imbalance in cell distribution between groups
    if max(group_counts) / sum(group_counts) > max_imbalance_ratio:
        return f"Skipped: imbalance between groups for AAV analysis"

    # Use raw counts
    subset_data.X = subset_data.layers['counts'].copy()

    # Filter genes expressed in at least 10 cells
    sc.pp.filter_genes(subset_data, min_cells=10)
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
## lowered min_cells_per_group to 30

def run_de_interaction_effects_analysis(adata, cell_type, min_cells_per_group=30, max_imbalance_ratio=0.9):
    # Filter data for the specific cell type
    subset_data = adata[adata.obs['neuron_cluster'] == cell_type].copy()

    # Checking for minimum cell count per group within the specific cell type
    group_counts = subset_data.obs.groupby(['AAV', 'TBI']).size()
    if any(group_counts < min_cells_per_group):
        return f"Skipped: not enough cells in one or more groups for interaction analysis"

    # Checking for imbalance in cell distribution between groups
    if max(group_counts) / sum(group_counts) > max_imbalance_ratio:
        return f"Skipped: imbalance between groups for interaction analysis"

    # Use raw counts
    subset_data.X = subset_data.layers['counts'].copy()

    # Filter genes expressed in at least 10 cells
    sc.pp.filter_genes(subset_data, min_cells=10)
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
    df_result['neuron_cluster'] = cell_type
    df_result['comparison'] = "TBI:AAV"

    return df_result