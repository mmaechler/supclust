\name{print.wilma}
\alias{print.wilma}
\title{Print Method for Wilma Objects}
\description{Yields an overview about the size and the final criterion
  values of the clusters that were selected by \code{wilma}.}
\usage{
print.wilma(x, ...)
}
\arguments{
  \item{x}{An \R object of \code{\link{class}} \code{"wilma"},
    typically the result of \code{\link{wilma}()}.}
  \item{\dots}{Further arguments passed to and from methods.}
}
\references{
  Marcel Dettling (2002)
  \emph{Supervised Clustering of Genes}, see
  \url{http://stat.ethz.ch/~dettling/supercluster.html}

  Marcel Dettling and Peter Bühlmann (2002).
  Supervised Clustering of Genes.
  \emph{Genome Biology}, \bold{3}(12): research0069.1-0069.15.
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{wilma}}}
\examples{
 ## Running the examples of Wilma's help page
 example(wilma, echo = FALSE)
 print(fit)
}
\keyword{classif}
\keyword{cluster}
