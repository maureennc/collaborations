# Code below is from 04-pairwise-select-volcano.R

library(readr)
library(EnhancedVolcano)
library(ggplot2)


##############################################################################

# IMPORT DATA

## Define paths
data_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/submission_01/tables/cell_type/'

data = read_csv(file.path(data_dir, "wald-test-cell_type-filtered.csv"))

##############################################################################

# CREATE DATA SUBSETS FOR FIGURE

AC <- data[data$comparison == "AC",]
CE <- data[data$comparison == "CE",]
DE <- data[data$comparison == "DE",]

##############################################################################

# AC - INHIBITORY NEURON

## Filter for "Excitatory neuron" cell type
sub_data <- AC[AC$cell_type == "Inhibitory neuron", ]
sub_data <- sub_data[sub_data$mean >= 0.1, ]

## Filter based on qval, log2fc, and mean expression
filtered_data <- sub_data[sub_data$qval < 0.05, ]
filtered_data <- filtered_data[abs(filtered_data$log2fc) >= 0.5, ]

## Extract filtered gene names
filtered_genes <- c('Meis1', 'Glra3', 'Rgs9', 'Rmst', 'Ptprd', 'Hs3st4', 'Unc13c', 'Sv2b',
                    'Gm48749', 'Dram1', 'Ring1', 'Ttr', 'Penk', 'Arhgap10', 'Aifm3')

## Volcano plot
AC_I <- EnhancedVolcano(sub_data,
                        lab = sub_data$gene,
                        x = 'log2fc',
                        y = 'pval',
                        pCutoffCol = 'qval',
                        title = paste("Inhibitory neuron - AC", sep=" - "),
                        pCutoff = 0.05,
                        FCcutoff = 0.5,
                        pointSize = 3.0,
                        labSize = 7,
                        selectLab = filtered_genes,
                        drawConnectors = TRUE,
                        widthConnectors = 0.75,
                        legendLabSize = 10,
                        max.overlaps = 20)

AC_I


##############################################################################

# AC - EXCITATORY NEURON

## Filter for "Excitatory neuron" cell type
sub_data <- AC[AC$cell_type == "Excitatory neuron", ]
sub_data <- sub_data[sub_data$mean >= 0.1, ]

## Filter based on qval, log2fc, and mean expression
filtered_data <- sub_data[sub_data$qval < 0.05, ]
filtered_data <- filtered_data[abs(filtered_data$log2fc) >= 0.5, ]

## Extract filtered gene names
filtered_genes <- 
  
  
  ## Volcano plot
  AC_E <- EnhancedVolcano(sub_data,
                          lab = sub_data$gene,
                          x = 'log2fc',
                          y = 'pval',
                          pCutoffCol = 'qval',
                          title = paste("Excitatory neuron - AC", sep=" - "),
                          pCutoff = 0.05,
                          FCcutoff = 0.5,
                          pointSize = 3.0,
                          labSize = 7,
                          selectLab = filtered_genes,
                          drawConnectors = TRUE,
                          widthConnectors = 0.75,
                          legendLabSize = 10,
                          max.overlaps = 20)

AC_E

##############################################################################

# CE - INHIBITORY NEURON

## Filter for "Excitatory neuron" cell type
sub_data <- CE[CE$cell_type == "Inhibitory neuron", ]
sub_data <- sub_data[sub_data$mean >= 0.1, ]

## Filter based on qval, log2fc, and mean expression
filtered_data <- sub_data[sub_data$qval < 0.05, ]
filtered_data <- filtered_data[abs(filtered_data$log2fc) >= 0.5, ]

## Extract filtered gene names
filtered_genes <- c('Rmst', 'Glra3', 'Tshz2', 'Sv2b', 'Pde7b', 'Ptprd', 'Nrg1', 'Penk', 'Cobl', 'Rgs9', 'Unc13c',
                    'Ttr', 'Prkcd', 'Zic4', 'Enpp2', 'Bcan')

## Volcano plot
CE_I <- EnhancedVolcano(sub_data,
                        lab = sub_data$gene,
                        x = 'log2fc',
                        y = 'pval',
                        pCutoffCol = 'qval',
                        title = paste("Inhibitory neuron - CE", sep=" - "),
                        pCutoff = 0.05,
                        FCcutoff = 0.5,
                        pointSize = 3.0,
                        labSize = 7,
                        selectLab = filtered_genes,
                        drawConnectors = TRUE,
                        widthConnectors = 0.75,
                        legendLabSize = 10,
                        max.overlaps = 20)

CE_I


##############################################################################

# CE - EXCITATORY NEURON

## Filter for "Excitatory neuron" cell type
sub_data <- CE[CE$cell_type == "Excitatory neuron", ]
sub_data <- sub_data[sub_data$mean >= 0.1, ]

## Filter based on qval, log2fc, and mean expression
filtered_data <- sub_data[sub_data$qval < 0.05, ]
filtered_data <- filtered_data[abs(filtered_data$log2fc) >= 0.5, ]

## Extract filtered gene names
filtered_genes <- c('Foxp2', 'Rorb', 'Prox1', 'Dpp10', 'Brinp3', 'Prox1',
                    'Ptprd', 'Disc1', 'Sv2c', 'Synpr', 'Rfx3', 'C1ql3', 'Rmst', 'Cobl', 'Gm18870')


## Volcano plot
CE_E <- EnhancedVolcano(sub_data,
                        lab = sub_data$gene,
                        x = 'log2fc',
                        y = 'pval',
                        pCutoffCol = 'qval',
                        title = paste("Excitatory neuron - CE", sep=" - "),
                        pCutoff = 0.05,
                        FCcutoff = 0.5,
                        pointSize = 3.0,
                        labSize = 7,
                        selectLab = filtered_genes,
                        drawConnectors = TRUE,
                        widthConnectors = 0.75,
                        legendLabSize = 10,
                        
                        max.overlaps = Inf,
                        xlim = c(-2,2))

CE_E

##############################################################################

# SAVE PLOTS

volcano_plots <- list(
  # "AC_inhibitory" = AC_I,
  #  "AC_excitatory" = AC_E,
  "CE_inhibitory" = CE_I,
  # "CE_excitatory" = CE_E
)

# save to new location in case the other ones change
#save_path <- "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/figures/wald/volcano-cell_type-curated"

for (i in seq_along(volcano_plots)) {
  plot_name <- paste0("volcano_", names(volcano_plots)[i], ".pdf")  # Use list name for the file
  ggsave(filename = file.path(save_path, plot_name), 
         plot = volcano_plots[[i]], 
         device = "pdf",  # Save as PDF
         width = 6, height = 6, dpi = 300)  # Adjust size and resolution as needed
}
