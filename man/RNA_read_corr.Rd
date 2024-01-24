\name{RNA_read_corr}
\alias{RNA_read_corr}
\title{Function to produce a scatterplot correlating RNA reads}
\description{
  This function produces a scatterplot correlating the number of detected genes with the total number of molecules identified in the transcriptome.
}
\usage{
  \code{RNA_read_corr(pbmc, cutoff_lo = 200, cutoff_hi = 2500)}
}
\arguments{
  \item{pbmc}{A Seurat object.}
  \item{cutoff_lo}{Lower cutoff for good-quality cells.}
  \item{cutoff_hi}{Upper cutoff for good-quality cells.}
}
\examples{
  \dontrun{
    scatterplot <- RNA_read_corr(pbmc, cutoff_lo = 200, cutoff_hi = 2500)
  }
}
\author{
  Your Name
}
\keyword{plot}
\export{RNA_read_corr}
