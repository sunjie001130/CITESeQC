\name{RNA_ADT_UMAP_corr}
\alias{RNA_ADT_UMAP_corr}
\title{Function to produce UMAP plots and scatterplot for RNA and ADT abundance}
\description{
  This function produces pairs of UMAP plots for the abundance of the input ADT and the corresponding gene expression. Additionally, it generates a scatterplot for the abundance of ADTs and the expression of the RNAs of the input gene.
}
\usage{
  \code{RNA_ADT_UMAP_corr(pbmc, rna_feature, adt_feature)}
}
\arguments{
  \item{pbmc}{A Seurat object.}
  \item{rna_feature}{Name of the RNA feature for plotting.}
  \item{adt_feature}{Name of the ADT feature for plotting.}
}
\examples{
  \dontrun{
    plots <- RNA_ADT_UMAP_corr(pbmc, rna_feature = "MYC", adt_feature = "CD45")
  }
}
\author{
  Your Name
}
\keyword{plot}
\export{RNA_ADT_UMAP_corr}
