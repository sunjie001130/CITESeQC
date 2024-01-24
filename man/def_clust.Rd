\name{def_clust}
\alias{def_clust}
\title{Function to define cell clusters using Seurat}
\description{
  This function defines cell clusters in a Seurat object using principal component analysis (PCA) and uniform manifold approximation and projection (UMAP).
}
\usage{
  \code{def_clust(pbmc, use_existing_clusters = FALSE, resolution = 0.8)}
}
\arguments{
  \item{pbmc}{A Seurat object.}
  \item{use_existing_clusters}{Logical, whether to use existing clusters.}
  \item{resolution}{Clustering resolution parameter.}
}
\examples{
  \dontrun{
    result <- def_clust(pbmc, use_existing_clusters = FALSE, resolution = 1.0)
  }
}
\author{
  Jie Sun
}
\keyword{clustering}
\keyword{export}{def_clust}
