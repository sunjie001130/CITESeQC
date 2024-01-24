\name{RNA_dist}
\alias{RNA_dist}
\title{Function to visualize specificity of gene expression across cell clusters}
\description{
  This function visualizes the specificity of the input gene expression across the cell clusters defined or imported through \code{def_clust()}.
}
\usage{
  \code{RNA_dist(pbmc, gene_name)}
}
\arguments{
  \item{pbmc}{A Seurat object.}
  \item{gene_name}{Name of the gene for visualization.}
}
\examples{
  \dontrun{
    plot <- RNA_dist(pbmc, gene_name = "MYC")
  }
}
\author{
  Your Name
}
\keyword{visualization}
\export{RNA_dist}
