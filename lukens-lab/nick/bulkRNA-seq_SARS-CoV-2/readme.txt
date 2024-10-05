Experimental design:
- Variable 1: Age (Young vs. Old)
- Variable 2: Infection (Mock-infected vs. SARS-CoV-2-infection)
- n = 4 mice per group; 16 mice total
- Sequencing: paired-end

Workflow:
- Step 1: Trimmed TruSeq2-PE adapters using Trimmomatic
- Step 2: Ran FastQC and MultiQC on samples
- Step 3: Used salmon to align reads to mouse genome (mm10)
- Step 4: Performed differential expression analysis in R using tximport and DESEq2

Sample information:
|Label_ID|Sample_ID|Age|Infection|
|T1|YM_1|Young|Mock|
|T2|YM_2|Young|Mock|
|T3|YM_3|Young|Mock|
|M4-Young|YM_4|Young|Mock|
|T4|YS_1|Young|SARS-CoV-2|
|T6|YS_2|Young|SARS-CoV-2|
|Y-MA10-33|YS_3|Young|SARS-CoV-2|
|Y-MA10-44|YS_4|Young|SARS-CoV-2|
|2Y-M1|AM_1|Aged|Mock|
|2Y-M2|AM_2|Aged|Mock|
|2Y-M3|AM_3|Aged|Mock|
|2Y-M4|AM_4|Aged|Mock|
|T8|AS_1|Aged|SARS-CoV-2|
|2Y-I1|AS_2|Aged|SARS-CoV-2|
|2Y-I2|AS_3|Aged|SARS-CoV-2|
|2Y-I3|AS_4|Aged|SARS-CoV-2|

Sample prefix identifiers:
- YM = young, mock
- YS = young, SARS-CoV-2-infected
- AM = aged, mock
- AS = aged, SARS-CoV-2-infected


Software version information:
- Trimmomatic 0.39
- FastQC 0.11.5
- MultiQC 1.11
- Salmon 1.5.1
- R 4.1.1

Output files:
- Salmon quant files
- FPKM matrix
- DESeq2 results tables (untransformed and apeglm-transformed for Log2FC visualization)
- DESeq2 normalized reads


Description of directories
- analysis: DESeq2 results output and annotated .Rmd analysis
- qc: FastQC and multiQC files
- raw_data: raw fastq files
- salmon_files: aligned files located in "quant_files" folder. Other folders are output for each salmon run
- salmon_mm10_index: contains mm10 reference genome for read alignment
- slurm_scripts: contains bash scripts for trimmomatic, fastqc, and salmon
- trimmed_paired: reads paired after trimming with trimmomatic
- trimmed_unpaired: reads discarded after trimming with trimmomatic

UVA Research Computing Acknowledgment statement:
-  Please cite Rivanna
- E.g. “The authors acknowledge Research Computing at The University of Virginia for providing computational resources and technical support that have contributed to the results reported within this publication. URL: https://rc.virginia.edu”
