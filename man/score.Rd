\name{score}
\alias{score}
\title{Wilcoxon Score for Binary Problems}
\description{For a set of \eqn{n} observations grouped into two classes (for
  example \eqn{n} expression values of a gene), the \code{score}
  function measures the separation of the classes. It can be interpreted
  as counting for each observation having response zero, the number of
  individuals of response class one that are smaller, and summing up
  these quantities.}
\usage{
score(x, resp)
}
\arguments{
  \item{x}{Numeric vector of length \eqn{n}, for example containing gene
    or cluster expression values of \eqn{n} different cases.}
  \item{resp}{Numeric vector of length \eqn{n} containing the \dQuote{binary}
    class labels of the cases.  Must be coded by \code{0} and \code{1}.}
}
\value{A numeric value, the \code{score}. The minimal \code{score} is
  zero, the maximal \code{score} is the product of the number of samples
  in class 0 and class 1. Values near the minimal or maximal
  \code{score} indicate good separation, whereas intermediate
  \code{score} means poor separation.}
\references{
  Marcel Dettling (2002)
  \emph{Supervised Clustering of Genes}, see
  \url{http://stat.ethz.ch/~dettling/supercluster.html}

  Marcel Dettling and Peter Bühlmann (2002).
  Supervised Clustering of Genes.
  \emph{Genome Biology}, \bold{3}(12): research0069.1-0069.15.
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{wilma}}, \code{\link{margin}} is the second statistic
  that is used there.
}
\examples{
data(leukemia, package="supclust")
op <- par(mfrow=c(1,3))
plot(leukemia.x[,69],leukemia.y)
title(paste("Score = ", score(leukemia.x[,69], leukemia.y)))

## Sign-flipping is very important
plot(leukemia.x[,161],leukemia.y)
title(paste("Score = ", score(leukemia.x[,161], leukemia.y),2))
x <- sign.flip(leukemia.x, leukemia.y)$flipped.matrix
plot(x[,161],leukemia.y)
title(paste("Score = ", score(x[,161], leukemia.y),2))
par(op)
}
\keyword{htest}
