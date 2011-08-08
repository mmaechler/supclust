\name{margin}
\alias{margin}
\title{Classification Margin Between Two Sample Classes}
\description{
  For a set of \eqn{n} observations grouped into two classes (for
  example \eqn{n} expression values of a gene), the \code{margin}
  function measures the size of the gap between the classes. This is the
  distance between the observation of response class zero having the
  lowest value, and the individual of with response one having the
  highest value.
}
\usage{
margin(x, resp)
}
\arguments{
  \item{x}{Numeric vector of length \eqn{n}, for example containing gene
    or cluster expression values of \eqn{n} different cases.}
  \item{resp}{Numeric vector of length \eqn{n} containing the \dQuote{binary}
    class labels of the cases.  Must be coded by \code{0} and \code{1}.}
}
\value{A numeric value, the \code{margin}. Positive \code{margin}
  indicates perfect separation of the response classes, whereas negative
  \code{margin} means imperfect separation.}
\references{
  Marcel Dettling (2002)
  \emph{Supervised Clustering of Genes}, see
  \url{http://stat.ethz.ch/~dettling/supercluster.html}

  Marcel Dettling and Peter Bühlmann (2002).
  Supervised Clustering of Genes.
  \emph{Genome Biology}, \bold{3}(12): research0069.1-0069.15.
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{wilma}}, \code{\link{score}} is the second statistic
  that is used there.
}
\examples{
data(leukemia, package="supclust")
op <- par(mfrow=c(1,3))
plot(leukemia.x[,69],leukemia.y)
title(paste("Margin = ", round(margin(leukemia.x[,69], leukemia.y),2)))

## Sign-flipping is very important
plot(leukemia.x[,161],leukemia.y)
title(paste("Margin = ", round(margin(leukemia.x[,161], leukemia.y),2)))
x <- sign.flip(leukemia.x, leukemia.y)$flipped.matrix
plot(x[,161],leukemia.y)
title(paste("Margin = ", round(margin(x[,161], leukemia.y),2)))
par(op)
}
\keyword{htest}
