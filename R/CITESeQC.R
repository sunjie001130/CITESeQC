#' Read CITE-seq data from an HDF5 file and create a Seurat object.
#'
#' @param filePath Path to the HDF5 file.
#' @return A Seurat object.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' }
#'
#' @import Seurat
#' @import hdf5r
#' @export
readCITEseqData <- function(filePath) {
  E = Read10X_h5(filename = filePath)
  pbmc <- CreateSeuratObject(counts = E$`Gene Expression`, project = "pbmc")
  pbmc[["ADT"]] <- CreateAssayObject(counts = E$`Antibody Capture`[, colnames(E$`Antibody Capture`) %in% colnames(pbmc)])
  return(pbmc)
}

#' Produce a scatterplot correlating the number of detected genes with the total number of molecules in the transcriptome.
#'
#' @param pbmc A Seurat object.
#' @param cutoff_lo Lower cutoff for good-quality cells.
#' @param cutoff_hi Upper cutoff for good-quality cells.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' RNA_read_corr(pbmc)
#' }
#'
#' @import Seurat
#' @export
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

#' Produce a scatterplot correlating the number of detected ADTs with the total number of ADT molecules on cell surfaces.
#'
#' @param pbmc A Seurat object.
#' @param cutoff_lo Lower cutoff for good-quality cells.
#' @param cutoff_hi Upper cutoff for good-quality cells.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' ADT_read_corr(pbmc)
#' }
#'
#' @import Seurat
#' @export
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

#' Produce a scatterplot correlating the number of molecules identified in the transcriptome with the percentage of mitochondrial genes.
#'
#' @param pbmc A Seurat object.
#' @param cutoff_lo Lower cutoff for good-quality cells.
#' @param cutoff_hi Upper cutoff for good-quality cells.
#' @param cutoff_mt Cutoff for mitochondrial percentage.
#'
#' @return A ggplot object and the updated Seurat object.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' pbmc <- RNA_mt_read_corr(pbmc)
#' }
#'
#' @import Seurat
#' @export
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

#' Filter cells based on specified cutoffs.
#'
#' @param pbmc A Seurat object.
#' @param feature_cutoff_lo Lower cutoff for the number of detected genes.
#' @param feature_cutoff_hi Upper cutoff for the number of detected genes.
#' @param mt_cutoff Cutoff for mitochondrial percentage.
#'
#' @return A trimmed Seurat object.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' pbmc.trim <- trimReads(pbmc)
#' }
#'
#' @import Seurat
#' @export
trimReads <- function(pbmc, feature_cutoff_lo = 200, feature_cutoff_hi = 2500, mt_cutoff = 10) {
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc.trim <- subset(pbmc, subset = nFeature_RNA > feature_cutoff_lo & nFeature_RNA < feature_cutoff_hi & percent.mt < mt_cutoff)
  return(pbmc.trim)
}

#' Define cell clusters using Seurat, either based on input gene expression matrix or by importing the definition.
#'
#' @param pbmc A Seurat object.
#' @param use_existing_clusters Logical, indicating whether to use existing clusters.
#' @param resolution Clustering resolution.
#'
#' @return A list containing the Seurat object and marker genes.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' result <- def_clust(pbmc)
#' }
#'
#' @import Seurat
#' @export
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

#' Visualize the specificity of input gene expression across cell clusters.
#'
#' @param pbmc A Seurat object.
#' @param gene_name Name of the gene for visualization.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' RNA_dist(pbmc, "your_gene_name")
#' }
#'
#' @import Seurat
#' @export
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

#' Histogram of normalized Shannon entropy values of marker genes identified in def_clust().
#'
#' @param pbmc A Seurat object.
#' @param markerListGEX List of marker genes.
#' @param top_n Number of top genes to consider.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' multiRNA_hist(pbmc, markerListGEX, top_n = 10)
#' }
#'
#' @import Seurat
#' @export
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

#' Visualize the specificity of input ADT abundance across cell clusters.
#'
#' @param pbmc A Seurat object.
#' @param gene_name Name of the ADT for visualization.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' ADT_dist(pbmc, "your_ADT_name")
#' }
#'
#' @import Seurat
#' @export
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

#' Histogram of normalized Shannon entropy values of all ADTs identified for cell clusters.
#'
#' @param pbmc A Seurat object.
#' @param ADTs List of ADTs.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' multiADT_hist(pbmc, ADTs)
#' }
#'
#' @import Seurat
#' @export
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

#' Scatterplot showing the correlation between the number of assayed genes in the transcriptome and the number of assayed cell surface proteins.
#'
#' @param pbmc A Seurat object.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' RNA_ADT_read_corr(pbmc)
#' }
#'
#' @import Seurat
#' @export
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

#' Pairs of UMAP plots and a scatterplot for the abundance of input ADT and the corresponding gene expression.
#'
#' @param pbmc A Seurat object.
#' @param rna_feature Name of the gene for gene expression.
#' @param adt_feature Name of the ADT for ADT abundance.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' RNA_ADT_UMAP_corr(pbmc, "your_gene_name", "your_ADT_name")
#' }
#'
#' @import Seurat
#' @export
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

#' Scatterplots for each cell cluster, showing the correlation between input ADT abundance and corresponding gene expression.
#'
#' @param pbmc A Seurat object.
#' @param rna_feature Name of the gene for gene expression.
#' @param adt_feature Name of the ADT for ADT abundance.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' RNA_ADT_cluster_corr(pbmc, "your_gene_name", "your_ADT_name")
#' }
#'
#' @import Seurat
#' @export
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

#' Histogram of the correlation coefficients in all pairs of ADTs and the corresponding genes in expression.
#'
#' @param seurat_object A Seurat object.
#' @param RNA_features Names of genes for gene expression.
#' @param ADT_features Names of ADTs for ADT abundance.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' RNA_ADT_hist(pbmc, c("gene1", "gene2"), c("ADT1", "ADT2"))
#' }
#'
#' @import Seurat
#' @export
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

#' Histograms showing the distribution of correlation coefficients in all pairs of ADTs and the corresponding genes for each cell cluster.
#'
#' @param seurat_object A Seurat object.
#' @param RNA_features Names of genes for gene expression.
#' @param ADT_features Names of ADTs for ADT abundance.
#'
#' @examples
#' \dontrun{
#' pbmc <- readCITEseqData("path/to/your/file.h5")
#' RNA_ADT_cluster_hist(pbmc, c("gene1", "gene2"), c("ADT1", "ADT2"))
#' }
#'
#' @import Seurat
#' @export
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







