\name{trimReads}
\alias{trimReads}
\title{Function to filter cells based on specified cutoffs}
\description{
  This function filters cells in a Seurat object based on user-defined cutoffs for the number of detected genes and mitochondrial percentage.
}
\usage{
  \code{trimReads(pbmc, feature_cutoff_lo = 200, feature_cutoff_hi = 2500, mt_cutoff = 10)}
}
\arguments{
  \item{pbmc}{A Seurat object.}
  \item{feature_cutoff_lo}{Lower cutoff for the number of detected genes.}
  \item{feature_cutoff_hi}{Upper cutoff for the number of detected genes.}
  \item{mt_cutoff}{Mitochondrial percentage cutoff.}
}
\examples{
  \dontrun{
    trimmed_pbmc <- trimReads(pbmc, feature_cutoff_lo = 300, feature_cutoff_hi = 2000, mt_cutoff = 5)
  }
}
\author{
  Your Name
}
\keyword{preprocessing}
\export{trimReads}
