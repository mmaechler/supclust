\name{leukemia}
\docType{data}
\alias{leukemia}
\alias{leukemia.x}
\alias{leukemia.y}
\alias{leukemia.z}
\title{A part of the Golub's famous AML/ALL-leukemia dataset}
\description{Part of the training set of the famous AML/ALL-leukemia
  dataset from the Whitehead Institute.  It has been reduced to 250
  genes, about half of which are very informative for classification,
  whereas the other half was chosen randomly.
}
\usage{
data(leukemia)
}
\format{Contains three \R-objects:
  The expression matrix \code{leukemia.x},
  the associated binary response variable \code{leukemia.y},
  and the associated 3-class response variable \code{leukemia.z}.
}
\source{\url{http://www.genome.wi.mit.edu/MPR/}}
\references{
  First published in\cr
  Golub et al. (1999)
  Molecular Classification of Cancer: Class Discovery and Class
  Prediction by Gene Expression Monitoring.
  \emph{Science} \bold{286}, 531--538.
}
\author{Marcel Dettling}
\examples{
data(leukemia, package="supclust")
str(leukemia.x)
str(leukemia.y)
str(leukemia.z)
op <- par(mfrow= 1:2)
plot(leukemia.x[,56], leukemia.y)
plot(leukemia.x[,174],leukemia.z)
par(op)
}
\keyword{datasets}



