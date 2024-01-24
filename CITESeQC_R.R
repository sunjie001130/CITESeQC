library(Seurat)
library(hdf5r)
library(tidyverse)
library(SeuratDisk)
library(patchwork)
library(grid)
library(gridExtra)
library(magick)
library(ggplot2)
library(celldex)
library(gridExtra)
# Read in data from the h5 file and create a Seurat object
readCITEseqData <- function(filePath) {
  E = Read10X_h5(filename = filePath)
  pbmc <- CreateSeuratObject(counts = E$`Gene Expression`, project = "pbmc")
  pbmc[["ADT"]] <- CreateAssayObject(counts = E$`Antibody Capture`[, colnames(E$`Antibody Capture`) %in% colnames(pbmc)])
  return(pbmc)
}

# ‘RNA_read_corr()’ produces a scatterplot correlating the number of detected genes with the total number of molecules identified in the transcriptome. 
# Since the cutoffs for good-quality cells will be passed as the arguments to the function, users can modify them for their data.  
RNA_read_corr <- function(pbmc, cutoff_lo = 200, cutoff_hi = 2500) {
  # Set suggested cutoffs
  feature_RNA_cutoff_lo <- cutoff_lo
  feature_RNA_cutoff_hi <- cutoff_hi
  data <- data.frame(
    x = pbmc$nFeature_RNA,
    y = pbmc$nCount_RNA
  )
  correlation <- cor.test(data$x, data$y)
  ggplot(data, aes(x = x, y = y)) +
    geom_point() +
    geom_vline(xintercept = c(feature_RNA_cutoff_lo, feature_RNA_cutoff_hi), color = "red", linetype = "dashed") +
    annotate("text", x = Inf, y = Inf, label = paste("Correlation:", round(correlation$estimate, 2), "; P value:", round(correlation$p.value, 2)), hjust = 1, vjust = 1) +
    xlim(c(min(data$x), max(data$x) + 1)) +
    labs(x = "nFeature_RNA", y = "nCount_RNA", title = "Scatter Plot with Correlation and Red Lines") +
    ggtitle("Correlation between the number of detected genes and the total number of molecules identified in the transcriptome") +
    theme_minimal()
}

# ‘ADT_read_corr()’ produces a scatterplot correlating the number of detected ADTs with the total number of ADT molecules identified on the cell surfaces. 
# Since the cutoffs identifying good-quality cells are annotated on the plot as passed as the arguments of the function, users can modify them for their data. 
ADT_read_corr <- function(pbmc, cutoff_lo = NULL, cutoff_hi = NULL) {
  # Set suggested cutoffs
  feature_RNA_cutoff_lo <- cutoff_lo
  feature_RNA_cutoff_hi <- cutoff_hi
  data <- data.frame(
    x = pbmc$nFeature_ADT,
    y = pbmc$nCount_ADT
  )
  correlation <- cor.test(data$x, data$y)
  ggplot(data, aes(x = x, y = y)) +
    geom_point() +
    geom_vline(xintercept = c(feature_RNA_cutoff_lo, feature_RNA_cutoff_hi), color = "red", linetype = "dashed") +
    annotate("text", x = Inf, y = Inf, label = paste("Correlation:", round(correlation$estimate, 2), "; P value:", round(correlation$p.value, 2)), hjust = 1, vjust = 1) +
    xlim(c(min(data$x), max(data$x) + 1)) +
    labs(x = "nFeature_ADT", y = "nCount_ADT", title = "Scatter Plot with Correlation and Red Lines") +
    ggtitle("Correlation between the number of detected ADTs and the total number of ADT molecules identified on the cell surfaces") +
    theme_minimal()
}

# ‘RNA_mt_read_corr()’ produces a scatterplot correlating the number of molecules identified in the transcriptome with the percentage of the mitochondrial genes. 
# With the correlation coefficient, users can test the hypothesis that the mitochondrial percentage should be constant regardless of the number of identified molecules.
RNA_mt_read_corr <- function(pbmc, cutoff_lo = 200, cutoff_hi = 2500, cutoff_mt = 10) {
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  # Set suggested cutoffs
  feature_RNA_cutoff_lo <- cutoff_lo
  feature_RNA_cutoff_hi <- cutoff_hi
  precent_MT_cutoff <- cutoff_mt
  data <- data.frame(
    x = pbmc$nFeature_RNA,
    y = pbmc$percent.mt
  )
  correlation <- cor.test(data$x, data$y)
  ggplot(data, aes(x = x, y = y)) +
    geom_point() +
    geom_vline(xintercept = c(cutoff_lo, cutoff_hi), color = "red", linetype = "dashed") +
    geom_hline(yintercept = precent_MT_cutoff, color = "red", linetype = "dashed") +
    annotate("text", x = Inf, y = Inf, label = paste("Correlation:", round(correlation$estimate, 2), "; P value:", round(correlation$p.value, 2)), hjust = 1, vjust = 1) +
    xlim(c(min(data$x), max(data$x) + 1)) +
    labs(x = "nFeature_RNA", y = "% of MT", title = "Scatter Plot with Correlation and Red Lines") +
    ggtitle("Correlation between the number of detected genes and the percentage of MT") +
    theme_minimal()
  return(pbmc)
}

# filter cells based the cutoff as input by user
trimReads <- function(pbmc, feature_cutoff_lo = 200, feature_cutoff_hi = 2500, mt_cutoff = 10) {
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc.trim <- subset(pbmc, subset = nFeature_RNA > feature_cutoff_lo & nFeature_RNA < feature_cutoff_hi & percent.mt < mt_cutoff)
  return(pbmc.trim)
}

# ‘def_clust()’ either defines the cell clusters based on the input gene expression matrix or imports the definition. 
# Specifically, it will employ Seurat to cluster the cells to the input clustering resolution and identify marker genes for each cluster for later use. 
def_clust <- function(pbmc, use_existing_clusters = FALSE, resolution = 0.8) {
    # Define new clusters
    DefaultAssay(pbmc) <- "RNA"
    pbmc <- NormalizeData(pbmc)
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc, features = all.genes)
    pbmc <- RunPCA(pbmc, npcs = 15, verbose = FALSE, features = VariableFeatures(object = pbmc))
    pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:15)
    pbmc <- FindNeighbors(pbmc, dims = 1:15, verbose = FALSE) 
    pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = resolution)
    Idents(pbmc) <- "seurat_clusters"
    
    # Generate a marker heatmap based upon the GEX clusters
    markerListGEX <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    markerListGEX %>%
      group_by(cluster) %>%
      top_n(n = 5, wt = avg_log2FC) -> top5
    DimPlot(pbmc, raster = FALSE, label = TRUE, label.size = 6, pt.size = 0.5) +
      ggtitle(sprintf("Clustering GEX Only | %d cells", ncol(pbmc))) + NoLegend()
    DoHeatmap(subset(pbmc, downsample = ceiling(0.02 * ncol(pbmc))), features = top5$gene) + NoLegend() +
      theme(axis.text = element_text(size = 5)) + ggtitle("GEX Markers in GEX Clustering")
    return(list(pbmc = pbmc, markerListGEX = markerListGEX))
}

# ‘RNA_dist()’ visualizes the specificity of the input gene expression across the cell clusters defined or 
# imported through def_clust(). For quantification and comparison, it will calculate Shannon entropy 
# on the expression distribution across the clusters 
RNA_dist <- function(pbmc, gene_name) {
  # Calculate gene expression for the specified gene
  gene_expression <- as.data.frame(AverageExpression(object = pbmc, assays = "RNA", features = gene_name, group.by = "ident", slot = "counts"))
  # Calculate relative frequencies for each cluster
  relative_frequencies <- t(gene_expression) / rowSums(gene_expression)
  # Calculate Shannon entropy for each cluster
  entropy <- -colSums(relative_frequencies * log2(relative_frequencies + 1e-10))  # Adding a small value to avoid log(0)
  # Plot
  plot(c(0:(ncol(gene_expression)-1)), t(gene_expression), type = "b", frame = FALSE, pch = 19, 
       col = "red", xlab = "Cluster", ylab = paste("Gene:", gene_name), ylim = c(0, 1.5)) 
  text(7, 1.5, paste("Shannon entropy for gene", gene_name, "across clusters:", round(entropy, 2)), pos = 1)
}


# ‘multiRNA_hist()’ is a histogram of normalized Shannon entropy values of the marker genes identified in def_clust().
# The histogram will plot how specific the marker genes are across the clusters. Users can modify the number of marker
# genes.  
multiRNA_hist <- function(pbmc, markerListGEX, top_n = 10) {
  # Get the top genes based on average log2FC
  top <- markerListGEX %>%
    group_by(cluster) %>%
    top_n(n = top_n, wt = avg_log2FC)
  # Initialize an empty vector to store the results
  results <- numeric(length(top$gene))
  # Loop through the list and perform calculations
  for (i in 1:length(top$gene)) {
    # Perform the calculations here
    gene_name <- top$gene[i] # Gene input by user
    gene_expression <- as.data.frame(AverageExpression(object = pbmc, assays = "RNA", features = gene_name, group.by = "ident", slot = "counts"))
    # Calculate relative frequencies for each cluster
    relative_frequencies <- t(gene_expression) / rowSums(gene_expression)
    # Calculate Shannon entropy for each cluster
    entropy <- -colSums(relative_frequencies * log2(relative_frequencies + 1e-10)) # Adding a small value to avoid log(0)
    # Store the result in the vector
    results[i] <- entropy
  }
  # Plot the histogram
  hist(results, 
       xlab = "Entropy",   # Specify the label for the x-axis
       main = "Histogram of Entropy of marker genes")
}

# ‘ADT_dist()’ visualizes the specificity of the input ADT abundance across the cell clusters. Specifically, it will calculate 
# Shannon entropy on the expression distribution across the clusters
ADT_dist <- function(pbmc, gene_name) {
  # Gene input by user
  gene_expression <- as.data.frame(AverageExpression(object = pbmc, assays = "ADT", features = gene_name, group.by = "ident", slot = "counts"))
  # Calculate relative frequencies for each cluster
  relative_frequencies <- t(gene_expression) / rowSums(gene_expression)
  # Calculate Shannon entropy for each cluster
  entropy <- -colSums(relative_frequencies * log2(relative_frequencies + 1e-10)) # Adding a small value to avoid log(0)
  # Plot
  plot(c(0:(ncol(gene_expression) - 1)), t(gene_expression), type = "b", frame = FALSE, pch = 19,
       col = "red", xlab = "Cluster", ylab = paste("ADT:", gene_name), ylim = c(0, 1600))
  text(7, 1550, paste("Shannon entropy for ADT", gene_name, "across clusters:", round(entropy, 2)), pos = 1)
}

# ‘multiADT_hist()’ is a histogram of normalized Shannon entropy values of all ADTs identified for the cell clusters. 
multiADT_hist <- function(pbmc, ADTs) {
  # Initialize an empty vector to store the results
  results <- numeric(length(ADTs))
  # Loop through the list and perform calculations
  for (i in 1:length(ADTs)) {
    # Perform the calculations here
    gene_name <- ADTs[i] # ADT input by user
    gene_expression <- as.data.frame(AverageExpression(object = pbmc, assays = "ADT", features = gene_name, group.by = "ident", slot = "counts"))
    # Calculate relative frequencies for each cluster
    relative_frequencies <- t(gene_expression) / rowSums(gene_expression)
    # Calculate Shannon entropy for each cluster
    entropy <- -colSums(relative_frequencies * log2(relative_frequencies + 1e-10)) # Adding a small value to avoid log(0)
    # Store the result in the vector
    results[i] <- entropy
  }
  # Plot histogram
  hist(results, 
       xlab = "Entropy",   # Specify the label for the x-axis
       main = "Histogram of Entropy of all ADTs in this dataset")
}

# ‘RNA_ADT_read_corr()’ produces a scatterplot showing the correlation between the number of assayed genes in the transcriptome and the number of assayed 
# cell surface proteins across the cells.
RNA_ADT_read_corr <- function(pbmc) {
  # Create a data frame
  data <- data.frame(
    x = pbmc$nCount_RNA,
    y = pbmc$nCount_ADT
  )
  # Calculate correlation
  correlation_count <- cor.test(data$x, data$y)
  # Plot correlation scatter plot
  ggplot(data, aes(x = x, y = y)) +
    geom_point() +  # Add scatter plot
    annotate("text", x = Inf, y = Inf, label = paste("Correlation:", round(correlation_count$estimate, 2), "; P value:", round(correlation_count$p.value, 2)), hjust = 1, vjust = 1) + # Add correlation text
    labs(x = "nCount_RNA", y = "nCount_ADT") + # Labels and title
    ggtitle("Correlation between the number of molecules in the transcriptome and in the cell surface proteome") +
    theme_minimal()  # Use a minimal theme
}

# ‘RNA_ADT_UMAP_corr()’ produces pairs of UMAP plots and a scatterplot. Each UMAP plot pair is drawn for the abundance of the input ADT and the corresponding gene expression, respectively. The scatterplot plots 
# the abundance of ADTs and the expression of the RNAs of the input gene.  
RNA_ADT_UMAP_corr <- function(pbmc, rna_feature, adt_feature) {
  # Feature plots for matched genes
  FeaturePlot(object = pbmc, features = rna_feature)
  FeaturePlot(object = pbmc, features = adt_feature)
  # calculate coefficient and p-value
  data <- FetchData(pbmc, vars = c(rna_feature, adt_feature))
  correlation <- cor.test(x = data[, 1], y = data[, 2], method = "spearman")
  # Plot feature scatter plot
  FeatureScatter(
    object = pbmc,
    feature1 = rna_feature,
    feature2 = adt_feature,
    cols = rep("black", 15)
  ) + NoLegend() +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste(
        "Correlation:",
        round(correlation$estimate, 2),
        "; P value:",
        round(correlation$p.value, 2)
      ),
      hjust = 1,
      vjust = 1
    ) +
    ggtitle("")
}

# ‘RNA_ADT_cluster_corr()’ is a set of scatterplots each drawn for each cell 
# cluster, showing the correlation between the input ADT abundance and the corresponding gene expression 
# for the cluster. 
RNA_ADT_cluster_corr <- function(pbmc, rna_feature, adt_feature) {
  Idents(pbmc) <- "seurat_clusters"
  # Get unique cluster IDs
  unique_clusters <- unique(pbmc$seurat_clusters)
  # Iterate through clusters and create feature plots
  plist <- list()
  for (cluster_id in levels(unique_clusters)) {
    # Subset the data for the current cluster
    subset_obj <- subset(pbmc, ident = cluster_id)
    # Calculate correlation
    data <- FetchData(subset_obj, vars = c(rna_feature, adt_feature))
    correlation <- cor.test(x = data[, 1], y = data[, 2], method = "spearman")
    # Create a feature plot for each input ADT
    feature_plot <- FeatureScatter(
      object = subset_obj,
      feature1 = rna_feature,
      feature2 = adt_feature,
      pt.size = 0.5
    ) +
      annotate(
        "text", 
        x = Inf,
        y = Inf,
        label = paste(
          "Correlation:",
          round(correlation$estimate, 2),
          "; P value:",
          round(correlation$p.value, 2)
        ),
        hjust = 1,
        vjust = 1
      ) +
      labs(title = paste("Cluster:", cluster_id), y = adt_feature) +
      NoLegend()
    
    plist[[cluster_id]] <- feature_plot
  }
  # Arrange the plots in a grid
  grid.arrange(grobs = plist, ncol = 4)
}

# ‘RNA_ADT_hist()’ is a histogram of the correlation coefficients in all pairs of ADTs 
# and the corresponding genes in expression.  
# Function to create scatter plots and histograms for correlation coefficients
RNA_ADT_hist <- function(seurat_object, RNA_features, ADT_features) {
  # Initialize a matrix to store correlation coefficients
  res_cor <- data.frame(matrix(ncol = 2, nrow = length(ADT_features)))
  names(res_cor) <- c("coef")
  # Loop through RNA and ADT features
  for (i in seq_along(ADT_features)) {
      # Fetch data and calculate correlation
      data <- FetchData(seurat_object, vars = c(RNA_features[i], ADT_features[i]))
      res_cor$coef[i] <- cor(x = data[, 1], y = data[, 2], method = "spearman")
  }
  # Plot histograms for correlation coefficients
  hist(res_cor$coef, 
       xlim = c(0, 1),
       xlab = "Coefficient",
       main = "Distribution of the correlation coefficients in pairs of ADTs and the corresponding genes")
}

# ‘RNA_ADT_cluster_hist()’ is a set of histograms, each showing the distribution of the correlation coefficients 
# in all pairs of ADTs and the corresponding genes for each cell cluster. 
RNA_ADT_cluster_hist <- function(seurat_object, RNA_features, ADT_features) {
  par(mfrow=c(4,4))
  par(mar=c(2,2,2,2))
  Idents(seurat_object) <- "seurat_clusters"
  unique_clusters <- unique(seurat_object$seurat_clusters)
  for (cluster_id in levels(unique_clusters)) {
    # Subset the data for the current cluster
    subset_obj <- subset(seurat_object, ident = cluster_id)
    # Initialize a matrix to store correlation coefficients
    res_cor <- data.frame(matrix(ncol = 1, nrow = length(ADT_features)))
    names(res_cor) <- c("coef")
    # Loop through RNA and ADT features
    for (i in seq_along(ADT_features)) {
      # Fetch data and calculate correlation
      data <- FetchData(seurat_object, vars = c(RNA_features[i], ADT_features[i]))
      res_cor$coef[i] <- cor(x = data[, 1], y = data[, 2], method = "spearman")
    }
    # Plot histograms for correlation coefficients
    hist(res_cor$coef, 
         xlab = "Coefficient", 
         main = paste("Cluster:", cluster_id))
  }
  par(mfrow=c(1,1))
}







