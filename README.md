# BF528_project5

This project aims to replicate the programmer and analyst task from project 4 using the data processed by data curator. The two major objectives of this replication study are to cluster the cells and identify the cell type using marker genes, then purpose novel gene markers for each cell type. 

Baron, Maayan, Adrian Veres, Samuel L. Wolock, Aubrey L. Faust, Renaud Gaujoux, Amedeo Vetere, Jennifer Hyoje Ryu, et al. 2016. “A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-Cell Population Structure.” Cell Systems 3 (4): 346–60.e4. PMID: 27667365

# Repository Contents
## Programmer.R
* Dependencies: R
* Inputs: quants_mat.gz
* Outputs: .rda file, plots for matrix features
* Tasks: Perform QC for UMI matrix, Identify clusters o cell type subpopulations

## Analyst.R
* Dependencies: R
* Inputs: .rda file
* Outputs: UMAP clustered image, heatmaps of marker expression, .csv of novel marker genes
* Tasks: Identify marker genes for each cluster, Label cluster cell types (manually), Visualize clustered cells, Visualize top marker genes per cluster, Identify novel marker genes
