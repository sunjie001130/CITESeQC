\name{multiRNA_hist}
\alias{multiRNA_hist}
\title{Function to create a histogram of normalized Shannon entropy values of marker genes}
\description{
  This function creates a histogram of normalized Shannon entropy values of the marker genes identified in \code{def_clust()}.
}
\usage{
  \code{multiRNA_hist(pbmc, markerListGEX, top_n = 10)}
}
\arguments{
  \item{pbmc}{A Seurat object.}
  \item{markerListGEX}{Marker gene list from \code{def_clust()}.}
  \item{top_n}{Number of top marker genes to consider.}
}
\examples{
  \dontrun{
    hist_plot <- multiRNA_hist(pbmc, markerListGEX, top_n = 15)
  }
}
\author{
  Jie Sun
}
\keyword{histogram}
\keyword{export}{multiRNA_hist}
