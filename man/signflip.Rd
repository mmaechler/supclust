\name{sign.flip}
\alias{sign.flip}
\title{Sign-flipping of Predictor Variables to Obtain Equal Polarity}
\description{Computes the \code{score} for each predictor variable
  (gene) in the \code{x}-Matrix, and multiplies its values with (-1) if
  its \code{score} is greater or equal than half of the maximal
  \code{score}. For gene expression data, this amounts to treating
  under- and overexpression symmetrically. After the sign-flip
  procedure, low (expression) values point towards response class 0 and
  high (expression) values point towards class 1.}
\usage{
sign.flip(x, y)
}
\arguments{
  \item{x}{Numeric matrix of explanatory variables (\eqn{p} variables in
    columns, \eqn{n} cases in rows). For example, these can be
    microarray gene expression data which should be sign-flipped and
    then clustered.}
  \item{y}{Numeric vector of length \eqn{n} containing the class labels
    of the individuals. These labels have to be coded by 0 and 1.}
}
\value{Returns a list containing:
  \item{flipped.matrix}{The sign-flipped \code{x}-matrix.}
  \item{signs}{Numeric vector of length \eqn{p}, which for each
    predictor variable indicates whether it was sign-flipped (coded by
    -1) or not (coded by +1).}
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
\seealso{\code{\link{wilma}} and \code{\link{score}}, as well as for a
  newer methodology, \code{\link{pelora}} and \code{\link{sign.change}}.}
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
par(op)# reset
}
\keyword{manip}
