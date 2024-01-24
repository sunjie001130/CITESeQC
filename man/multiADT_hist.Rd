\name{multiADT_hist}
\alias{multiADT_hist}
\title{Function to create a histogram of normalized Shannon entropy values of all ADTs}
\description{
  This function creates a histogram of normalized Shannon entropy values of all ADTs identified for the cell clusters.
}
\usage{
  \code{multiADT_hist(pbmc, ADTs)}
}
\arguments{
  \item{pbmc}{A Seurat object.}
  \item{ADTs}{Vector of ADTs for which entropy is calculated.}
}
\examples{
  \dontrun{
    hist_plot <- multiADT_hist(pbmc, ADTs = c("CD45", "CD3", "CD19"))
  }
}
\author{
  Your Name
}
\keyword{histogram}
\export{multiADT_hist}
