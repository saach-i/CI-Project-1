# De-convoluting cancer genome evolution following TP53 mutation in high grade serous ovarian cancer.
## CI Project 1

This repository contains the code required to reproduce the results presented in the following paper:

S.Sachdev, H.Lu, S.Ghaem-Maghami, De-convoluting cancer genome evolution following TP53 mutation in high grade serous ovarian cancer (2023).

# Data
All the data files required to reproduce the results in this paper are stored in the Data folder. This study used  publicly available datasets from the TCGA HGSOC study, previously published here https://www.nature.com/articles/nature10166.

TCGA Ploidy data contains the ploidy data, the output from ABSOLUTE.
TCGA_ov_clinical contains the clinical data for this patient cohort, including survival outcomes.
TCGA.OV.sampleMap_SNP6_nocnv_genomicSegment contains the segmented copy number data.
limma_signficantCN_TCGA contains the all-lesions output from the GISTIC 2.0 analysis.
mutations_ov_tcga contains the mutation data for the key HGSOC genes.
TCGA_OV_MCP contains the output from the MCPcounter estimation, including TLS levels.

# Code
All RScripts used in this study are stored in the Script folder.

