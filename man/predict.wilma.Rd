\name{predict.wilma}
\alias{predict.wilma}
\title{Predict Method for Wilma}
\description{Yields fitted values or predicted class labels for training
  and test data, which are based on the supervised gene clusters
  \code{wilma} found, and on a choice of four different classifiers: the
  nearest-neighbor rule, diagonal linear discriminant analysis, logistic
  regression and aggregated trees.}
\usage{
predict.wilma(object, newdata = NULL, type = c("fitted", "class"),
              classifier = c("nnr", "dlda", "logreg", "aggtrees"),
              noc = object$noc, ...)
}
\arguments{
  \item{object}{An \R object of \code{\link{class}} \code{"wilma"},
    typically the result of \code{\link{wilma}()}.}
  \item{newdata}{Numeric matrix with the same number of explanatory
    variables as the original \code{x}-matrix (\eqn{p} variables in
    columns, \eqn{r} cases in rows). For example, these can be
    additional microarray gene expression data which should be
    predicted.}
  \item{type}{Character string, describing whether fitted values
    \code{"fitted"} or predicted class labels \code{"class"} should be
    returned.}
  \item{classifier}{Character string, saying which classifier should be
    used. Choices are \code{"nnr"}, the 1-nearest-neighbor-rule;
    \code{"dlda"}, diagonal linear discriminant analysis;
    \code{"logreg"}, logistic regression; \code{"aggtrees"} aggregated
    trees. }
  \item{noc}{Integer, saying with how many clusters the fitted values or
    class label predictions should be determined. Also numeric vectors
    are allowed as an argument. The output is then a numeric matrix with
    fitted values or class label predictions for a multiple number of
    clusters.} 
  \item{\dots}{}
}
\details{
  If \code{newdata = NULL}, then the in-sample fitted values or class
  label predictions are returned. 
}
\value{
  Depending on whether \code{noc} is a single number or a numeric
  vector. In the first case, a numeric vector of length \eqn{r} is
  returned, which contains fitted values for \code{noc} clusters, or
  class label predictions with \code{noc} clusters.

  In the latter case, a numeric matrix with \code{length(noc)} columns,
  each containing fitted values for \code{noc} clusters, or class label
  predictions with \code{noc} clusters, is returned.
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
\seealso{\code{\link{wilma}} and for the four classifiers,
  \code{\link{nnr}}, \code{\link{dlda}}, \code{\link{logreg}},
  \code{\link{aggtrees}}.}
\examples{
## Working with a "real" microarray dataset
data(leukemia, package="supclust")

## Generating random test data: 3 observations and 250 variables (genes)
set.seed(724)
xN <- matrix(rnorm(750), nrow = 3, ncol = 250)

## Fitting Wilma
fit  <- wilma(leukemia.x, leukemia.y, noc = 3, trace = 1)

## Fitted values and class predictions for the training data
predict(fit, type = "cla")
predict(fit, type = "fitt")

## Predicting fitted values and class labels for test data
predict(fit, newdata = xN)
predict(fit, newdata = xN, type = "cla", classifier = "nnr", noc = c(1,2,3))
predict(fit, newdata = xN, type = "cla", classifier = "dlda", noc = c(1,3))
predict(fit, newdata = xN, type = "cla", classifier = "logreg")
predict(fit, newdata = xN, type = "cla", classifier = "aggtrees")
}
\keyword{classif}
\keyword{cluster}
