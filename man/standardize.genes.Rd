\name{standardize.genes}
\alias{standardize.genes}
\title{Standardization of Predictor Variables}
\description{Standardizes each column (gene) of the \code{x}-matrix to
  zero mean and unit variance. This function is not to be called by the
  user, the standardization is handled internally in \code{pelora}.}
\usage{
standardize.genes(exmat)
}
\arguments{
  \item{exmat}{Numeric matrix of explanatory variables (\eqn{p}
    variables in columns, \eqn{n} cases in rows). For example, these can
    be microarray gene expression data which should be standardized and
    then grouped.}
}
\value{
  Returns a list containing:
  \item{x}{The standardized \code{x}-matrix}
  \item{means}{Numeric vector of length \eqn{p}, containing the initial
    mean of each column (gene) of the \code{x}-matrix.}
  \item{sdevs}{Numeric vector of length \eqn{p}, containing the initial
    standard deviation of each column (gene) of the \code{x}-matrix.}
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
\seealso{\code{\link{pelora}}}
\keyword{manip}
