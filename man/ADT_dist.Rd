\name{ADT_dist}
\alias{ADT_dist}
\title{Function to visualize specificity of ADT abundance across cell clusters}
\description{
  This function visualizes the specificity of the input ADT abundance across the cell clusters defined or imported through \code{def_clust()}.
}
\usage{
  \code{ADT_dist(pbmc, gene_name)}
}
\arguments{
  \item{pbmc}{A Seurat object.}
  \item{gene_name}{Name of the ADT for visualization.}
}
\examples{
  \dontrun{
    plot <- ADT_dist(pbmc, gene_name = "CD45")
  }
}
\author{
  Jie Sun
}
\keyword{visualization}
\keyword{export}
