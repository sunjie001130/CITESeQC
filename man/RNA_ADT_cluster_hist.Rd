\name{RNA_ADT_cluster_hist}
\alias{RNA_ADT_cluster_hist}
\title{Function to create histograms for correlation coefficients of RNA and ADT pairs in each cell cluster}
\description{
  This function is a set of histograms, each showing the distribution of the correlation coefficients in all pairs of ADTs and the corresponding genes for each cell cluster.
}
\usage{
  \code{RNA_ADT_cluster_hist(seurat_object, RNA_features, ADT_features)}
}
\arguments{
  \item{seurat_object}{A Seurat object.}
  \item{RNA_features}{Vector of RNA features for correlation.}
  \item{ADT_features}{Vector of ADT features for correlation.}
}
\examples{
  \dontrun{
    hist_plots <- RNA_ADT_cluster_hist(pbmc, RNA_features = c("MYC", "CD45"), ADT_features = c("CD3", "CD19"))
  }
}
\author{
  Jie Sun
}
\keyword{histogram}
\export{RNA_ADT_cluster_hist}
