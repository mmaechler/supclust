\name{coef.pelora}
\alias{coef.pelora}
\title{Extract the Model Coefficients of Pelora}
\description{Yields the coefficients of the penalized logistic
  regression model that is fitted by \code{pelora} with its groups of
  predictor variables (genes) as input}
\usage{
coef.pelora(object, ...)
}
\arguments{
  \item{object}{An \R object of \code{\link{class}} \code{"pelora"},
    typically the result of \code{\link{pelora}()}.}
  \item{\dots}{Further arguments passed to and from methods.}
}
\value{A numeric vector of length \eqn{\code{noc}+1}, giving the
  penalized logistic regression coefficients for the intercept and the
  \code{noc} groups and/or single variables identified by \code{pelora}.}
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
 coef(fit)
}
\keyword{classif}
\keyword{cluster}
