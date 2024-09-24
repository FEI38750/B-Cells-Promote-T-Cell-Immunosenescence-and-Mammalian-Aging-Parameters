# Code deposition for the paper: B Cells Promote T Cell Immunosenescence and Mammalian Aging Parameters
1. *CellRanger_commands:* Commands used for preprocessing single-cell transcriptomics fastq files.
2. *B_T_scRNAseq_analysis.R:* This R script performs 5’ single-cell transcriptomics, including VDJ T cell receptor analyses on splenic CD3+ cells from aged µMT and WT mice. The gene lists used for volcano plot are in the folder 'volcano_plot_genes'.
3. *PICseq_10XSingleCell.R:* This R script is used for analyzing single-cell RNA-seq data from B and T cell and physically interacting cells (PICs).
4. *B2T_nichenet.R and T2B_nichenet.R:* Those scripts utilizes the NicheNet framework to analyze B-to-T cell interactions, identifying ligand-receptor pairs that may mediate communication between these cells.
5. *PICseq_script.R:* This R script is specifically used for PICseq analysis.
6. *Classifier4PICs.ipynb:* This Jupyter notebook is used for training and applying a Stacked Ensemble model to predict PICs. The matrices required for training and prediction are imported from the script PICseq_script.R between lines 84 and 99.
7. *GeoMx_code_Spleen.rtf:* Code for processing GeoMx data.

