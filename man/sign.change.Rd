\name{sign.change}
\alias{sign.change}
\title{Sign-flipping of Predictor Variables to Obtain Equal Polarity}
\description{Computes the empirical correlation for each predictor
  variable (gene) in the \code{x}-Matrix with the response \code{y}, and
  multiplies its values with (-1) if the empirical correlation has a
  negative sign. For gene expression data, this amounts to treating
  under- and overexpression symmetrically. After the \code{sign.change},
  low (expression) values point towards response class 0 and high
  (expression) values point towards class 1.}
\usage{
sign.change(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric matrix of explanatory variables (\eqn{p} variables in
    columns, \eqn{n} cases in rows). For example, these can be
    microarray gene expression data which should be sign-flipped and
    then grouped.}
  \item{y}{Numeric vector of length \eqn{n} containing the class labels
    of the individuals. These labels have to be coded by 0 and 1.}
}
\value{
  Returns a list containing:
  \item{x.new}{The sign-flipped \code{x}-matrix.}
  \item{signs}{Numeric vector of length \eqn{p}, which for each
    predictor variable indicates whether it was sign-flipped (coded by
    -1) or not (coded by +1).}
}
\references{
  Marcel Dettling (2003)
  \emph{Finding Groups of Predictive Genes from Microarray Data}, see
  \url{http://stat.ethz.ch/~dettling/supervised.html}

  Marcel Dettling and Peter Bühlmann (2004).
  Finding Predictive Gene Groups from Microarray Data.
  To appear in the \emph{Journal of Multivariate Analysis}.
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{pelora}}, as well as for older methodology,
  \code{\link{wilma}} and \code{\link{sign.flip}}.}
\examples{
data(leukemia, package="supclust")

op <- par(mfrow=c(1,3))
plot(leukemia.x[,69],leukemia.y)
title(paste("Margin = ", round(margin(leukemia.x[,69], leukemia.y),2)))

## Sign-flipping is very important
plot(leukemia.x[,161],leukemia.y)
title(paste("Margin = ", round(margin(leukemia.x[,161], leukemia.y),2)))
x <- sign.change(leukemia.x, leukemia.y)$x.new
plot(x[,161],leukemia.y)
title(paste("Margin = ", round(margin(x[,161], leukemia.y),2)))
par(op)
}
\keyword{manip}
