\name{readCITEseqData}
\alias{readCITEseqData}
\title{Function to read CITE-seq data}
\description{
  This function reads CITE-seq data from a file and returns a Seurat object.
}
\usage{
  \code{readCITEseqData(filePath)}
}
\arguments{
  \item{filePath}{Path to the CITE-seq data file.}
}
\value{
  A Seurat object.
}
\examples{
  \dontrun{
    pbmc <- readCITEseqData("path/to/your/file.h5")
  }
}
\author{
  Jie Sun
}
\keyword{data}
\keyword{export}
