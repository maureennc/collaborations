PLOS Pathogens submission:
Caspase-1 in Cx3cr1-expressing cells drives an IL-18-dependent T cell response that promotes parasite control during acute Toxoplasma gondii infection


Link to T. Gondii-infected spleen dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207173

The publicly available count matrix from GSE207173 (Clark et al., 2023) was imported into a Python environment and processed using Scanpy. For quality thresholds, cells with less than 1k counts, cells greater than 25k counts, or cells containing greater than 20% mitochondrial RNA content were filtered and excluded from the analysis. The data matrix was normalized and logarithmized, and PCA was used for dimensionality reduction using the ARPACK wrapper. Unsupervised clustering was performed using the leiden algorithm, with a resolution set to 0.6. Differential expression was performed using Wilcoxon rank sum test.

Scripts

- 01_GSE207173: Performs QC filtering, library size normalization, PCA, dimensionality reduction, and manual cell type annotation based on leiden cluster differential expression results.

- 02_GSE207173: Performs differential expression of CX3CR1- and CX3CR1+ cells to address reviewer comment: "Characterize the population of Caspase-1 expressing Cx3cr1+ cells."

- Performed GO overrepresentation analysis using DAVID

- 03_GSE207173: Used DAVID results in bubble chart data visualization.