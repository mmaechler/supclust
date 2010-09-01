\name{coef.pelora}
\alias{coef.pelora}
\title{Extract the Model Coefficients of Pelora}
\description{Yields the coefficients of the penalized logistic
  regression model that is fitted by \code{pelora} with its groups of
  predictor variables (genes) as input}
\usage{
\S3method{coef}{pelora}(object, ...)
}
\arguments{
  \item{object}{an \R object of \code{\link{class}} \code{"pelora"},
    typically the result of \code{\link{pelora}()}.}
  \item{\dots}{further arguments passed to and from methods.}
}
\value{A numeric vector of length \eqn{\code{noc}+1}, giving the
  penalized logistic regression coefficients for the intercept and the
  \code{noc} groups and/or single variables identified by \code{pelora}.}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{pelora}}, also for references.}
\examples{
 ## Running the examples of Pelora's help page
 example(pelora, echo = FALSE)
 coef(fit)
}
\keyword{classif}
\keyword{cluster}
