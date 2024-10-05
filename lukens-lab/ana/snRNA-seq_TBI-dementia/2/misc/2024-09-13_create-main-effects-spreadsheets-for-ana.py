#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 10:46:59 2024

@author: maureen
"""

# CREATE TABLES FOR MAIN EFFECTS FOR NEURONS

import os
import pandas as pd

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/multivariate"

###############################################################################

# Read csv

filtered_aav = pd.read_csv(os.path.join(save_dir, 'aav-main-effects-cell_type-filtered.csv'))

filtered_tbi = pd.read_csv(os.path.join(save_dir, 'tbi-main-effects-cell_type-filtered.csv'))

###############################################################################

# Filter

groups = ['Excitatory neuron', 'Inhibitory neuron']

filtered_aav = filtered_aav[filtered_aav['cell_type'].isin(groups)]
filtered_tbi = filtered_tbi[filtered_tbi['cell_type'].isin(groups)]


## Significant genes only
aav_sig = filtered_aav[filtered_aav['qval'] < 0.05].copy()
tbi_sig = filtered_tbi[filtered_tbi['qval'] < 0.05].copy()


## Top hits
aav_top = aav_sig[(aav_sig['mean'] > 0.1) & (abs(aav_sig['log2fc']) > 0.5)].copy()
tbi_top = tbi_sig[(tbi_sig['mean'] > 0.1) & (abs(tbi_sig['log2fc']) > 0.5)].copy()



###############################################################################

# Export top hits as csv

csv_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/de/multivariate/top-hits"

aav_top.to_csv(os.path.join(csv_dir, 'aav-main-effects-top-hits.csv'), index = False)
tbi_top.to_csv(os.path.join(csv_dir, 'tbi-main-effects-top-hits.csv'), index = False)

filtered_aav.to_csv(os.path.join(csv_dir, 'aav-main-effects-results.csv'), index = False)
filtered_tbi.to_csv(os.path.join(csv_dir, 'tbi-main-effects-results.csv'), index = False)


###############################################################################