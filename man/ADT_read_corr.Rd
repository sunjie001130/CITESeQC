\name{ADT_read_corr}
\alias{ADT_read_corr}
\title{Function to produce a scatterplot correlating ADT reads}
\description{
  This function produces a scatterplot correlating the number of detected ADTs with the total number of ADT molecules identified on the cell surfaces.
}
\usage{
  \code{ADT_read_corr(pbmc, cutoff_lo = NULL, cutoff_hi = NULL)}
}
\arguments{
  \item{pbmc}{A Seurat object.}
  \item{cutoff_lo}{Lower cutoff for good-quality cells.}
  \item{cutoff_hi}{Upper cutoff for good-quality cells.}
}
\examples{
  \dontrun{
    scatterplot <- ADT_read_corr(pbmc, cutoff_lo = 100, cutoff_hi = 2000)
  }
}
\author{
  Your Name
}
\keyword{plot}
\export{ADT_read_corr}
