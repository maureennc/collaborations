# Set working directory

dir = '/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/global-sting-ko/results/spreadsheets/de'
setwd(dir)
getwd()

library(readr)

#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

##############################################################################

# IMPORT DATA

cell_type <- read_csv('global-sting-wt-ko-DE-cell_type.csv')

## Create data subsets
microglia <- cell_type[cell_type$cell_type == 'Microglia', ]
astrocyte <- cell_type[cell_type$cell_type == 'Astrocyte', ]
excitatory <- cell_type[cell_type$cell_type == 'Excitatory neuron', ]
inhibitory <- cell_type[cell_type$cell_type == 'Inhibitory neuron', ]
oligodendrocyte <- cell_type[cell_type$cell_type == 'Oligodendrocyte', ]
opc <- cell_type[cell_type$cell_type == 'OPC', ]
endothelial <- cell_type[cell_type$cell_type == 'Endothelial cell', ]
#pericyte <- cell_type[cell_type$cell_type == 'Pericyte', ]

microglia_filtered <- microglia[microglia$qval < 0.05 & abs(microglia$log2fc) >= 0.5, ]
astrocyte_filtered <- astrocyte[astrocyte$qval < 0.05 & abs(astrocyte$log2fc) >= 0.5, ]
excitatory_filtered <- excitatory[excitatory$qval < 0.05 & abs(excitatory$log2fc) >= 0.5, ]
inhibitory_filtered <- inhibitory[inhibitory$qval < 0.05 & abs(inhibitory$log2fc) >= 0.5, ]
oligodendrocyte_filtered <- oligodendrocyte[oligodendrocyte$qval < 0.05 & abs(oligodendrocyte$log2fc) >= 0.5, ]
opc_filtered <- opc[opc$qval < 0.05 & abs(opc$log2fc) >= 0.5, ]
endothelial_filtered <- endothelial[endothelial$qval < 0.05 & abs(endothelial$log2fc) >= 0.5, ]
#pericyte_filtered <- pericyte[pericyte$qval < 0.05 & abs(pericyte$log2fc) >= 0.5, ]

##############################################################################

# VOLCANO PLOT

## Microglia

microglia_genes = c('Prr5l', 'Sox5', 'Stat1', 'Tns1', 'Dpp10', 'Ifi204', 'Ptgds', 'Fcrls', 'Opcml', 'Axl', 'Apobec1', 'Cck', 'Cd83', 'Clu', 'Ext1', 'Cnksr2',
                    'Slc1a2', 'Camk2a', 'Nrgn', 'Brinp1')

p1 <- EnhancedVolcano(microglia,
                      lab = microglia$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'Microglia',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 6,
                      #selectLab = microglia_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE, 
                      widthConnectors = 0.5,
                      #xlim = c(-2, 2),
                      ylim = c(0, 8),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = 15
)
p1

##############################################################################

## Excitatory neuron

excitatory_genes = c('Ebf2', 'Lef1', 'Tcf7l2', 'Edaradd', 'Ebf1', 'Baiap3', 'Btbd11', 'Erbb4', 'Adarb2', 'Wnt3', 'Rasef', 'Spp1', 
                     'Igf1', 'Ak7', 'Agt', 'Nos1', 'Reln', 'Tspan2', 'Lgr5', 'Notch3',
                     'Myo5n', 'Ttr', 'Fap', 'Olfr111', 'Arc', 'Inf2',  'Vxn', 'Npas2', 'Jph1',
                     'Pparg', 'Etv5', 'Irf3', 'Il1rap',
                     'Ovol2', 'Satb2', 'Satb1', 'Epha4',
                     'Map3k5', 'Map3k1', 'Camk2a', 'Igf1')


p2 <- EnhancedVolcano(excitatory,
                      lab = excitatory$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'Excitatory neurons',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 6,
                      selectLab = excitatory_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      #xlim = c(-2, 2),
                      #ylim = c(0, 10),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = 60
)

p2

##############################################################################

## Inhibitory neuron

inhibitory_genes = c('Ttr', 'Cpa6', 'Kcnq5', 'Scg2', 'Celf2', 'Ano3', 'Lingo2', 'Foxp1', 'Sox5', 'Slit1', 'Inf2', 'Mef2c', 'Camk2a', 'Otx2os1', 'Lmo1',
                     'Dlk1', 'Slc17a6', 'Oxt', 'Nek10', 'Notch3', 'Col24a1', 'Acvr1c', 'Mybpc1', 'Abcb5', 'Jcad')

p3 <- EnhancedVolcano(inhibitory,
                      lab = inhibitory$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'Inhibitory neurons',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 6,
                      selectLab = inhibitory_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.50,
                      #xlim = c(-3, 3),
                      #ylim = c(0, 17),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      #legendLabSize = 10,
                      max.overlaps = 25
                      
)

p3

##############################################################################

#astrocyte_genes = c('Ttr')

p4 <- EnhancedVolcano(astrocyte,
                      lab = astrocyte$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'Astrocytes',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 6,
                      #selectLab = astrocyte_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      #xlim = c(-2, 2),
                      #ylim = c(0, 17),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = 20
                      
)

p4


##############################################################################

## Oligodendrocytes

oligo_genes = c('Rgs20', 'Gli2', 'Ryr3', 'Gm3764', 'Ndst3', 'Ak5', 'Sel1l3', 'Vgf', 'Hif3a', 'Ctsd', 'Nell1', 'Ttr', 'Alk',
                'Grn', 'Alk', 'Hs6st3', 'Lingo1', 'Sv2b', 'Slit3')

p5 <- EnhancedVolcano(oligodendrocyte,
                      lab = oligodendrocyte$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'Oligodendrocytes',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 6,
                      #selectLab = oligo_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      xlim = c(-2, 2),
                      #ylim = c(0, 17),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = 15
                      
)

p5

##############################################################################

## OPC

opc_genes = c('Cst3', 'Wwtr1', 'Ccdc114', 'A230006K03Rik', 'Ddx60', 'Smpd3', 'Hcn1', 'Airn')

p6 <- EnhancedVolcano(opc,
                      lab = opc$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'OPC',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 6,
                      #selectLab = opc_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      #xlim = c(-2, 2),
                      #ylim = c(0, 17),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = 15
                      
)

p6

##############################################################################

## Endothelial cells

p7 <- EnhancedVolcano(endothelial,
                      lab = endothelial$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'Endothelial cells',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 6,
                      #selectLab = endothelial_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      #xlim = c(-2, 2),
                      #ylim = c(0, 6),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = 20
                      
)

p7

##############################################################################

# Save the plots
setwd("~/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/global-sting-ko/results/figures/volcano")
plot_list <- list(p1, p2, p3, p4, p5, p6, p7)

plot_names <- c("Microglia", "Excitatory_neuron", "Inhibitory_neuron",
                "Astrocyte", "Oligodendrocyte", "OPC", "Endothelial")

## Save as pdf
for (i in seq_along(plot_list)) {
  ggsave(filename = paste0(plot_names[i], ".pdf"),
         plot = plot_list[[i]],
         width = 6, height = 6, dpi = 300,
         path = getwd())
}

## Save as png
for (i in seq_along(plot_list)) {
  ggsave(filename = paste0(plot_names[i], ".png"),
         plot = plot_list[[i]],
         width = 6, height = 6, dpi = 1000,
         path = getwd())
}


##############################################################################
