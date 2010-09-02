\name{print.wilma}
\alias{print.wilma}
\title{Print Method for Wilma Objects}
\description{Yields an overview about the size and the final criterion
  values of the clusters that were selected by \code{wilma}.}
\usage{
\S3method{print}{wilma}(x, ...)
}
\arguments{
  \item{x}{An \R object of \code{\link{class}} \code{"wilma"},
    typically the result of \code{\link{wilma}()}.}
  \item{\dots}{Further arguments passed to and from methods.}
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{wilma}}, also for references.}
\examples{
 ## Running the examples of Wilma's help page
 example(wilma, echo = FALSE)
 print(fit)
}
\keyword{classif}
\keyword{cluster}
