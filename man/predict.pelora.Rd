\name{predict.pelora}
\alias{predict.pelora}
\title{Predict Method for Pelora}
\description{Yields fitted values, predicted class labels and
  conditional probability estimates for training and test data, which
  are based on the gene groups \code{pelora} found, and on its internal
  penalized logistic regression classifier.}
\usage{
predict.pelora(object, newdata = NULL, newclin = NULL,
               type = c("fitted", "probs", "class"), noc = object$noc, ...)
}
\arguments{
  \item{object}{An \R object of \code{\link{class}} \code{"pelora"},
    typically the result of \code{\link{pelora}()}.}
  \item{newdata}{Numeric matrix with the same number of explanatory
    variables as the original \code{x}-matrix (\eqn{p} variables in
    columns, \eqn{r} cases in rows). For example, these can be
    additional microarray gene expression data which should be
    predicted.}  
  \item{newclin}{Numeric matrix with the same number of additional
    (clinical) explanatory variables as the original \code{u}-matrix
    (\eqn{m} variables in columns, \eqn{r} cases in rows) that are used 
    in the (penalized logistic regression) prediction model, but neither
    grouped nor averaged. Only needs to be given, if the model fit
    included an \code{u}-matrix. For example, these can be 'traditional'
    clinical variables.}
  \item{type}{Character string, describing whether fitted values
    \code{"fitted"}, estimated conditional probabilites \code{"probs"}
    or class labels \code{"class"} should be returned.}
  \item{noc}{Integer, saying with how many clusters the fitted values,
    probability estimates or class labels should be determined. Also
    numeric vectors are allowed as an argument. The output is then a
    numeric matrix with fitted values, probability estimates or class
    labels for a multiple number of clusters.}
  \item{\dots}{Further arguments passed to and from methods.}
}
\details{
  If \code{newdata = NULL}, then the in-sample fitted values,
  probability estimates and class label predictions are returned.
}
\value{
  Depending on whether \code{noc} is a single number or a numeric
  vector. In the first case, a numeric vector of length \eqn{r} is
  returned, which contains fitted values for \code{noc} clusters, or
  probability estimates/class label predictions with \code{noc}
  clusters.

  In the latter case, a numeric matrix with \code{length(noc)} columns,
  each containing fitted values for \code{noc} clusters, or
  probability estimates/class label predictions with \code{noc}
  clusters, is returned.}
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
## Working with a "real" microarray dataset
data(leukemia, package="supclust")

## Generating random test data: 3 observations and 250 variables (genes)
set.seed(724)
xN <- matrix(rnorm(750), nrow = 3, ncol = 250)

## Fitting Pelora
fit <- pelora(leukemia.x, leukemia.y, noc = 3)

## Fitted values and class probabilities for the training data
predict(fit, type = "cla")
predict(fit, type = "prob")

## Predicting fitted values and class labels for the random test data
predict(fit, newdata = xN)
predict(fit, newdata = xN, type = "cla", noc = c(1,2,3))
predict(fit, newdata = xN, type = "pro", noc = c(1,3))

## Fitting Pelora such that the first 70 variables (genes) are not grouped
fit <- pelora(leukemia.x[, -(1:70)], leukemia.y, leukemia.x[,1:70])

## Fitted values and class probabilities for the training data
predict(fit, type = "cla")
predict(fit, type = "prob")

## Predicting fitted values and class labels for the random test data
predict(fit, newdata = xN[, -(1:70)], newclin = xN[, 1:70])
predict(fit, newdata = xN[, -(1:70)], newclin = xN[, 1:70], "cla", noc  = 1:10)
predict(fit, newdata = xN[, -(1:70)], newclin = xN[, 1:70], type = "pro")
}
\keyword{classif}
\keyword{cluster}
