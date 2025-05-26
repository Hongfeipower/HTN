# Single Cell Multi-Omics Data Reveals Heterogeneity in Liver Tissue Microenvironment Induced by Hypertension

## Abstract
Hypertension (HTN) may induce liver damage, yet the effects on liver cell subpopulations remain obscure. Understanding these microenvironmental changes could offer early HTN and liver disease indicators. We employed single-cell multi-omics and histone chromatin immunoprecipitation sequencing (ChIP-seq) to scrutinize microenvironmental alterations between normal and HTN liver tissues. Analysis revealed an HTN-related hepatocyte subpopulation, termed “Hepatocytes_1”, via single-cell RNA sequencing (scRNA-seq). Five potential pathogenic genes (UPB1, SDS, PCCA, CYP3A4, and PPARGC1A) were identified in Hepatocytes_1. Additionally, the regulatory network of cis-regulatory elements (CREs) in genes within Hepatocytes_1 was unveiled using single-cell ATAC sequencing (scATAC-seq) and histone ChIP-seq. Specifically, enhancers in HTN cells recruited the HTN-related transcription factor (TF) MLXIPL. The active promoter of the disease-associated gene CYP3A4 formed a TF-mediated regulatory network, resulting in heightened expression in HTN cells. These findings of HTN-related gene markers and CREs in the liver offer novel insights for HTN prevention and drug targeting.
![Graphical Abstract](https://github.com/Hongfeipower/HTN/blob/main/Graphical%20Abstract.png)

## Folder Description
- ATAC
  * ATAC_GeneScore_AverageExpression.csv is a CSV file recording GeneScore for scATAC-seq.
  * peak.csv is a CSV file containing all peak positions and cell types.
- RNA
  * Hepatocytes_1_eG.csv
  * Hepatocytes_1_grup_differential_gene.csv
  * cell_type.csv
  * new_type_all_markers.csv
  * new_type_number_cluster_group.csv
  * original_all_markers_intergrate.csv
  * original_number_cluster_group.csv

- Code
  * Step 1. scRNA.r: R code for analyzing scATAC-seq data using ArchR.
  * Step 2. Pseudotime.R: R code for trajectory analysis.
  * Step 3. constrained_scATAC.R: R code for analyzing scATAC-seq data using ArchR.
  * Step 4. new_network.R :Code for constructing promoter and enhancer regulatory networks.
  
- promoter_enhancer:bedfile for promoter and enhancer
  * HTN_Hep1_enhancer.bed is an enhancer BED file for the HTN state of Hepatocytes_1.
  * HTN_Hep1_promoter.bed is an promoter BED file for the HTN state of Hepatocytes_1.
  * normal_Hep1_enhancer.bed is an enhancer BED file for the normal state of Hepatocytes_1.
  * normal_Hep1_promoter.bed is an promoter BED file for the normal state of Hepatocytes_1.
