\name{dlda}
\alias{dlda}
\alias{nnr}
\alias{logreg}
\alias{aggtrees}
\title{Classification with Wilma's Clusters}
\description{The four functions \code{nnr} (nearest neighbor rule),
  \code{dlda} (diagonal linear discriminant analysis), \code{logreg}
  (logistic regression) and \code{aggtrees} (aggregated trees) are used
  for binary classification with the cluster representatives of Wilma's
  output.
}
\usage{
dlda    (xlearn, xtest, ylearn)
nnr     (xlearn, xtest, ylearn)
logreg  (xlearn, xtest, ylearn)
aggtrees(xlearn, xtest, ylearn)
}

\arguments{
  \item{xlearn}{Numeric matrix of explanatory variables (\eqn{q}
    variables in columns, \eqn{n} cases in rows), containing the learning
    or training data. Typically, these are the (gene) cluster
    representatives of Wilma's output.}
  \item{xtest}{A numeric matrix of explanatory variables (\eqn{q}
    variables in columns, \eqn{m} cases in rows), containing the test or
    validation data. Typically, these are the fitted (gene) cluster
    representatives of Wilma's output for the training data, obtained
    from \code{predict.wilma}.}
  \item{ylearn}{Numeric vector of length \eqn{n} containing the class labels
    for the training observations. These labels have to be coded by 0 and 1.}
}
\details{\code{nnr} implements the 1-nearest-neighbor-rule with
  Euclidean distance function. \code{dlda} is linear discriminant
  analysis, using the restriction that the covariance matrix is diagonal
  with equal variance for all predictors. \code{logreg} is default
  logistic regression. \code{aggtrees} fits a default stump (a
  classification tree with two terminal nodes) by \code{rpart} for every
  predictor variable and uses majority voting to determine the final
  classifier.
}
\value{Numeric vector of length \eqn{m}, containing the predicted class
  labels for the test observations. The class labels are coded by 0 and
  1.}
\references{see those in \code{\link{wilma}}.}% maintain in *one* place
\author{Marcel Dettling}% no longer \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{wilma}}}
\examples{
## Generating random learning data: 20 observations and 10 variables (clusters)
set.seed(342)
xlearn <- matrix(rnorm(200), nrow = 20, ncol = 10)

## Generating random test data: 8 observations and 10 variables(clusters)
xtest  <- matrix(rnorm(80),  nrow = 8,  ncol = 10)

## Generating random class labels for the learning data
ylearn <- as.numeric(runif(20)>0.5)

## Predicting the class labels for the test data
nnr(xlearn, xtest, ylearn)
dlda(xlearn, xtest, ylearn)
logreg(xlearn, xtest, ylearn)
aggtrees(xlearn, xtest, ylearn)
}
\keyword{classif}
