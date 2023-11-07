# CITESeQC

CITESeQC is an R package that enables comprehensive, quantitative, and reproducible quality control for CITE-seq data.

1. For comprehensive quality control, CITESeQC assesses the quality of not only the gene expression layer and the surface protein layer but also their interactions.
2. For quantitative QC, CITESeQC quantifies the specificity of marker genes or protein markers across clusters and their correlational relationships using Spearman’s correlation and visualizes the relationship distribution across cell clusters in the histogram.
3. For reproducible QC, CITESeQC will guide users through a simple linear process that results in a full markdown report with supporting figures and explanations with minimal intervention from the user.

##  Installation
1. Install the package using the "devtools" package
   
```
# Install devtools package if not already installed
if (!require(devtools)) {
  install.packages("devtools")
}

# Install the package from GitHub
devtools::install_github("sunjie001130/CITESeQC")
```


2. Install the package using the "remotes" package
   
```
# Install remotes package if not already installed
if (!require(remotes)) {
  install.packages("remotes")
}

# Install your package from GitHub
remotes::install_github("sunjie001130/CITESeQC")
```

##  Tutorial
An example R markdown file in PDF version, as well as the dataset used in the example, are provided.

## Library Dependencies
Please install the following packages before using the CITESeQC tool.
- library(Seurat)
- library(hdf5r)
- library(tidyverse)
- library(SeuratDisk)
- library(patchwork)
- library(grid)
- library(gridExtra)
- library(magick)
- library(ggplot2)
- library(celldex)
