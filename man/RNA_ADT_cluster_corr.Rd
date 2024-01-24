\name{RNA_ADT_cluster_corr}
\alias{RNA_ADT_cluster_corr}
\title{Function to produce scatterplots for RNA and ADT abundance in each cell cluster}
\description{
  This function is a set of scatterplots, each drawn for each cell cluster, showing the correlation between the input ADT abundance and the corresponding gene expression for the cluster.
}
\usage{
  \code{RNA_ADT_cluster_corr(pbmc, rna_feature, adt_feature)}
}
\arguments{
  \item{pbmc}{A Seurat object.}
  \item{rna_feature}{Name of the RNA feature for plotting.}
  \item{adt_feature}{Name of the ADT feature for plotting.}
}
\examples{
  \dontrun{
    scatterplots <- RNA_ADT_cluster_corr(pbmc, rna_feature = "MYC", adt_feature = "CD45")
  }
}
\author{
  Jie Sun
}
\keyword{plot}
\keyword{export}{RNA_ADT_cluster_corr}
