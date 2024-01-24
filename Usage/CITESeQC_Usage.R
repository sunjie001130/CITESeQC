# How to use the CITESeQC tool

# Read in data from the h5 file and create a Seurat object
filePath <- "/Users/sunjie/Desktop/CiteSeq/CITESeQC/pbmc_10k_protein_v3_filtered_feature_bc_matrix.h5"
pbmc <- readCITEseqData(filePath)

# ‘RNA_read_corr()’
RNA_read_corr(pbmc)

# ‘ADT_read_corr()’
ADT_read_corr(pbmc)

# ‘RNA_mt_read_corr()’
RNA_mt_read_corr(pbmc)

# If user wants to filter out bad cells, run this step
pbmc.trim <- trimReads(pbmc)

# ‘def_clust()’
pbmc_clust <- def_clust(pbmc.trim)

# ‘RNA_dist()’
RNA_dist(pbmc_clust$pbmc, "CCR7")

# ‘multiRNA_hist()’
multiRNA_hist(pbmc_clust$pbmc, pbmc_clust$markerListGEX)

# ‘ADT_dist()’
ADT_dist(pbmc_clust$pbmc, "CD14-TotalSeqB")

# multiADT_hist()
ADTs <- c("CD14-TotalSeqB", "CD4-TotalSeqB", "TIGIT-TotalSeqB", "CD8a-TotalSeqB", "CD19-TotalSeqB")
multiADT_hist(pbmc_clust$pbmc, ADTs)

# ‘RNA_ADT_read_corr()’
RNA_ADT_read_corr(pbmc_clust$pbmc)

# ‘RNA_ADT_UMAP_corr()’
RNA_ADT_UMAP_corr(pbmc_clust$pbmc, "rna_CD14", "adt_CD14-TotalSeqB")

# RNA_ADT_cluster_corr()
RNA_ADT_cluster_corr(pbmc_clust$pbmc, "rna_CD14", "adt_CD14-TotalSeqB")

# RNA_ADT_hist()
RNA_features <- c('rna_CD14', 'rna_CD4', 'rna_TIGIT', 'rna_CD8A', 'rna_CD19')
ADT_features <- c('adt_CD14-TotalSeqB', 'adt_CD4-TotalSeqB', 'adt_TIGIT-TotalSeqB', 'adt_CD8a-TotalSeqB', 'adt_CD19-TotalSeqB')
RNA_ADT_hist(pbmc_clust$pbmc, RNA_features, ADT_features)


# RNA_ADT_cluster_hist()
RNA_ADT_cluster_hist(pbmc_clust$pbmc, RNA_features, ADT_features)


