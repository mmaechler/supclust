\name{print.pelora}
\alias{print.pelora}
\title{Print Method for Pelora Objects}
\description{Yields an overview about the type, size and final criterion
  value of the predictor variables that were selected by \code{pelora}.}
\usage{
print.pelora(x, digits = getOption("digits"), details = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An \R object of \code{\link{class}} \code{"pelora"},
    typically the result of \code{\link{pelora}()}.}
  \item{digits}{The number of digits that should be printed.}
  \item{details}{Logical, defaults to \code{FALSE}. If set to
    \code{TRUE}, the output corresponds to
    \code{\link{summary.pelora}}.}
  \item{\dots}{Further arguments passed to and from methods.}
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
\examples{
 ## Running the examples of Pelora's help page
 example(pelora, echo = FALSE)
 print(fit)
}
\keyword{classif}
\keyword{cluster}
