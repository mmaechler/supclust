\name{plot.pelora}
\alias{plot.pelora}
\title{2-Dimensional Visualization of Pelora's Output}
\description{Yields a projection of the cases (for example \eqn{n} gene
  expression profiles) into the space of the first two gene group
  centroids that were identified by \code{pelora}.}
\usage{
\S3method{plot}{pelora}(x, main = "2-Dimensional Projection Pelora's output",
            xlab = NULL, ylab = NULL, col = seq(x$yvals), ...)
}
\arguments{
  \item{x}{An \R object of \code{\link{class}} \code{"pelora"},
    typically the result of \code{\link{pelora}()}.}
  \item{main}{A character string, giving the title of the plot.}
  \item{xlab}{A character string, giving the annotation of the
    \code{x}-axis.}
  \item{ylab}{A character string, giving the annotation of the
    \code{x}-axis.}
  \item{col}{A numeric vector of length 2, coding the colors that will
    be used for plotting the class labels.}
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
 plot(fit)
}
\keyword{classif}
\keyword{cluster}

