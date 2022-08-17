# SPS_2022

### _A novel survival prediction signature outperforms PAM50 and artificial intelligence-based feature-selection methods_

The relevant raw datasets were obtained from the following sources:
1. RNAseq data "_CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz_" from Cancer Cell Line Encyclopedia (https://data.broadinstitute.org/ccle/)
2. BRCA survival metadata "_model_list_20220124.csv_" from Cell Model Passports (https://cellmodelpassports.sanger.ac.uk/downloads)
3. BRCA subtype metadata "_Expression_22Q1_Public_subsetted.csv_" from DepMap Portal (https://depmap.org/portal/download/custom/)
4. Three independent clinical datasets from GSE81538, GSE202203, and TCGA-BRCA available publicly online
 
---
 
Datasets 1-3 were first pre-processed in **ccle-preprocessing.R** before main processing in **ccle-processing-main.R**.

Dataset 4 was used for translating labels to clinical datasets via **label_transfer_1(Calculate_Dis).R** and **label_transfer_2(MarkLabel).ipynb**

Boruta codes and derived features can be found in "boruta" folder.

RGES results can be found in "gene_reversal_score_outputs" folder, using codes from Bin et al. (https://github.com/Bin-Chen-Lab/RGES).
