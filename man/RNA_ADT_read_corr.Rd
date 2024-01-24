\name{RNA_ADT_read_corr}
\alias{RNA_ADT_read_corr}
\title{Function to produce a scatterplot correlating RNA and ADT reads}
\description{
  This function produces a scatterplot showing the correlation between the number of assayed genes in the transcriptome and the number of assayed cell surface proteins across the cells.
}
\usage{
  \code{RNA_ADT_read_corr(pbmc)}
}
\arguments{
  \item{pbmc}{A Seurat object.}
}
\examples{
  \dontrun{
    scatterplot <- RNA_ADT_read_corr(pbmc)
  }
}
\author{
  Jie Sun
}
\keyword{plot}
\keyword{export}
