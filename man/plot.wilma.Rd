\name{plot.wilma}
\alias{plot.wilma}
\title{2-Dimensional Visualization of Wilma's Output}
\description{Yields a projection of the cases (for example \eqn{n} gene 
  expression profiles) into the space of the first two gene group
  centroids that were identified by \code{wilma}.}
\usage{
plot.wilma(x, xlab = NULL, ylab = NULL, col = seq(x$yvals),
           main = "2-Dimensional Projection of Wilma's Output", ...)
}
\arguments{
  \item{x}{An \R object of \code{\link{class}} \code{"wilma"},
    typically the result of \code{\link{wilma}()}.}
  \item{xlab}{A character string, giving the annotation of the
    \code{x}-axis.}
  \item{ylab}{A character string, giving the annotation of the
    \code{x}-axis.}
  \item{col}{A numeric vector of length 2, coding the colors that will
    be used for plotting the class labels.}
  \item{main}{A character string, giving the title of the plot.}
  \item{\dots}{Further arguments passed to and from methods.}
}
\references{
  Marcel Dettling (2002)
  \emph{Supervised Clustering of Genes}, see
  \url{http://stat.ethz.ch/~dettling/supercluster.html}

  Marcel Dettling and Peter Bühlmann (2002).
  Supervised Clustering of Genes.
  \emph{Genome Biology}, \bold{3}(12): research0069.1-0069.15.
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{wilma}}}
\examples{
 ## Running the examples of Wilma's help page
 example(wilma, echo = FALSE)
 plot(fit)
}
\keyword{classif}
\keyword{cluster}
