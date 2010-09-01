\name{print.pelora}
\alias{print.pelora}
\title{Print Method for Pelora Objects}
\description{Yields an overview about the type, size and final criterion
  value of the predictor variables that were selected by \code{pelora}.}
\usage{
\S3method{print}{pelora}(x, digits = getOption("digits"), details = FALSE, ...)
}
\arguments{
  \item{x}{an \R object of \code{\link{class}} \code{"pelora"},
    typically the result of \code{\link{pelora}()}.}
  \item{digits}{the number of digits that should be printed.}
  \item{details}{logical, defaults to \code{FALSE}. If set to
    \code{TRUE}, the output corresponds to
    \code{\link{summary.pelora}}.}
  \item{\dots}{Further arguments passed to and from methods.}
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{pelora}}, also for references.
}
\examples{
 ## Running the examples of Pelora's help page
 example(pelora, echo = FALSE)
 print(fit)
}
\keyword{classif}
\keyword{cluster}
