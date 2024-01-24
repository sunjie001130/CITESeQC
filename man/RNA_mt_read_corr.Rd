\name{RNA_mt_read_corr}
\alias{RNA_mt_read_corr}
\title{Function to produce a scatterplot correlating RNA reads and mitochondrial percentage}
\description{
  This function produces a scatterplot correlating the number of detected genes with the percentage of mitochondrial genes in the transcriptome.
}
\usage{
  \code{RNA_mt_read_corr(pbmc, cutoff_lo = 200, cutoff_hi = 2500, cutoff_mt = 10)}
}
\arguments{
  \item{pbmc}{A Seurat object.}
  \item{cutoff_lo}{Lower cutoff for good-quality cells.}
  \item{cutoff_hi}{Upper cutoff for good-quality cells.}
  \item{cutoff_mt}{Mitochondrial percentage cutoff.}
}
\examples{
  \dontrun{
    scatterplot <- RNA_mt_read_corr(pbmc, cutoff_lo = 200, cutoff_hi = 2500, cutoff_mt = 5)
  }
}
\author{
  Jie Sun
}
\keyword{plot}
\keyword{export}
