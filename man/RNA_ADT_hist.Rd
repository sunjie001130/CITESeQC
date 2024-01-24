\name{RNA_ADT_hist}
\alias{RNA_ADT_hist}
\title{Function to create histograms for correlation coefficients of RNA and ADT pairs}
\description{
  This function creates a histogram of the correlation coefficients in all pairs of ADTs and the corresponding genes in expression.
}
\usage{
  \code{RNA_ADT_hist(seurat_object, RNA_features, ADT_features)}
}
\arguments{
  \item{seurat_object}{A Seurat object.}
  \item{RNA_features}{Vector of RNA features for correlation.}
  \item{ADT_features}{Vector of ADT features for correlation.}
}
\examples{
  \dontrun{
    hist_plot <- RNA_ADT_hist(pbmc, RNA_features = c("MYC", "CD45"), ADT_features = c("CD3", "CD19"))
  }
}
\author{
  Jie Sun
}
\keyword{histogram}
\keyword{export}{RNA_ADT_hist}
