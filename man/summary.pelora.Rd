\name{summary.pelora}
\alias{summary.pelora}
\title{Summary Method for Pelora Objects}
\description{Yields detailed information about the variables (genes)
  that have been selected, and how they were grouped.}
\usage{
\S3method{summary}{pelora}(object, digits, \dots)
}
\arguments{
  \item{object}{an \R object of \code{\link{class}} \code{"pelora"},
    typically the result of \code{\link{pelora}()}.}
  \item{digits}{The number of digits that should be printed.}
  \item{\dots}{Further arguments passed to and from methods.}
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{pelora}}}
\examples{
 ## Running the examples of Pelora's help page
 example(pelora, echo = FALSE)
 summary(fit)
}
\keyword{classif}
\keyword{cluster}
