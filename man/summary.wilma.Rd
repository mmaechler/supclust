\name{summary.wilma}
\alias{summary.wilma}
\title{Summary Method for Wilma Objects}
\description{Yields detailed information about the variables (genes)
  that have been selected, and how they were clustered.}
\usage{
\S3method{summary}{wilma}(object, ...)
}
\arguments{
  \item{object}{An \R object of \code{\link{class}} \code{"wilma"},
    typically the result of \code{\link{wilma}()}.}
  \item{\dots}{Further arguments passed to and from methods.}
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{wilma}}, also for references.}
\examples{
 ## Running the examples of Wilma's help page
 example(wilma, echo = FALSE)
 summary(fit)
}
\keyword{classif}
\keyword{cluster}
