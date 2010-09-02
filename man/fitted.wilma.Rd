\name{fitted.wilma}
\alias{fitted.wilma}
\title{Extract the Fitted Values of Wilma}
\description{Yields the fitted values, i.e. the centroids of the (gene)
  clusters that have been found by \code{wilma}.}
\usage{
\S3method{fitted}{wilma}(object, \dots)
}
\arguments{
  \item{object}{An \R object of \code{\link{class}} \code{"wilma"},
    typically the result of \code{\link{wilma}()}.}
  \item{\dots}{further arguments passed to and from methods.}
}
\value{
  Numeric matrix of fitted values (for \eqn{n} cases in rows, and
  \code{noc} group centroids in columns).
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{wilma}}, also for references.}
\examples{
 ## Running the examples of Wilma's help page
 example(wilma, echo = FALSE)
 fitted(fit)
}
\keyword{classif}
\keyword{cluster}
