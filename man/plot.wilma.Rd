\name{plot.wilma}
\alias{plot.wilma}
\title{2-Dimensional Visualization of Wilma's Output}
\description{Yields a projection of the cases (for example \eqn{n} gene
  expression profiles) into the space of the first two gene group
  centroids that were identified by \code{wilma}.}
\usage{
\S3method{plot}{wilma}(x, xlab = NULL, ylab = NULL, col = seq(x$yvals),
           main = "2-Dimensional Projection of Wilma's Output", ...)
}
\arguments{
  \item{x}{an \R object of \code{\link{class}} \code{"wilma"},
    typically the result of \code{\link{wilma}()}.}
  \item{xlab}{character string, giving the annotation of the
    \code{x}-axis.}
  \item{ylab}{character string, giving the annotation of the
    \code{x}-axis.}
  \item{col}{a numeric vector of length 2, coding the colors that will
    be used for plotting the class labels.}
  \item{main}{a character string, giving the title of the plot.}
  \item{\dots}{Further arguments passed to and from methods.}
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{wilma}}, also for references.}
\examples{
 ## Running the examples of Wilma's help page
 example(wilma, echo = FALSE)
 plot(fit)
}
\keyword{classif}
\keyword{cluster}
