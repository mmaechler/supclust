\name{wilma}
\alias{wilma}
\title{Supervised Clustering of Predictor Variables}
\description{Performs supervised clustering of predictor variables for
  large (microarray gene expression) datasets. Works in a greedy forward
  strategy and optimizes a combination of the Wilcoxon and Margin
  statistics for finding the clusters.}
\usage{
wilma(x, y, noc, genes = NULL, flip = TRUE, once.per.clust = FALSE, trace = 0)
}
\arguments{
  \item{x}{Numeric matrix of explanatory variables (\eqn{p} variables in
    columns, \eqn{n} cases in rows). For example, these can be
    microarray gene expression data which should be clustered.}
  \item{y}{Numeric vector of length \eqn{n} containing the class labels
    of the individuals. These labels have to be coded by 0 and 1.}
  \item{noc}{Integer, the number of clusters that should be searched for
    on the data.} 
  \item{genes}{Defaults to \code{NULL}. An optional list (of length
    \code{noc}) of vectors containing the indices (column numbers) of
    the previously known initial clusters.} 
  \item{flip}{Logical, defaults to \code{TRUE}. Is indicating whether
    the clustering should be done with or without sign-flipping.}
  \item{once.per.clust}{Logical, defaults to \code{FALSE}. Is indicating
    if each variable (gene) should only be allowed to enter into each
    cluster once; equivalently, the cluster mean profile has only
    weights \eqn{\pm 1}{+/- 1} for each variable.}
  \item{trace}{Integer >= 0; when positive, the output of the internal
    loops is provided; \code{trace >= 2} provides output even from the
    internal C routines.}
}

\value{
  \code{wilma} returns an object of class "wilma". The functions
  \code{print} and \code{summary} are used to obtain an overview of the
  clusters that have been found. The function \code{plot} yields a
  two-dimensional projection into the space of the first two clusters
  that \code{wilma} found. The generic function \code{fitted} returns
  the fitted values, these are the cluster representatives. Finally,
  \code{predict} is used for classifying test data on the basis of
  Wilma's cluster with either the nearest-neighbor-rule, diagonal linear
  discriminant analysis, logistic regression or aggregated trees.

  An object of class "wilma" is a list containing:
  \item{clist}{A list of length \code{noc}, containing integer vectors
    consisting of the indices (column numbers) of the variables (genes)
    that have been clustered.} 
  \item{steps}{Numerical vector of length \code{noc}, showing the number
    of forward/backward cycles in the fitting process of each cluster.}
  \item{y}{Numeric vector of length \eqn{n} containing the class labels
    of the individuals. These labels have to be coded by 0 and 1.}
  \item{x.means}{A list of length \code{noc}, containing numerical
    matrices consisting of the cluster representatives after insertion
    of each variable.}
  \item{noc}{Integer, the number of clusters that has been searched for
    on the data.}
  \item{signs}{Numerical vector of length \eqn{p}, saying whether the
    \eqn{i}th variable (gene) should be sign-flipped (-1) or not (+1).}
}
\references{
  Marcel Dettling (2002)
  \emph{Supervised Clustering of Genes}, see
  \url{http://stat.ethz.ch/~dettling/supercluster.html}

  Marcel Dettling and Peter Bühlmann (2002).
  Supervised Clustering of Genes.
  \emph{Genome Biology}, \bold{3}(12): research0069.1-0069.15.

  Marcel Dettling and Peter Bühlmann (2004).
  Finding Predictive Gene Groups from Microarray Data.
  To appear in the \emph{Journal of Multivariate Analysis}.
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{score}}, \code{\link{margin}}, and for a newer
  methodology, \code{\link{pelora}}.}
\examples{
## Working with a "real" microarray dataset
data(leukemia, package="supclust")

## Generating random test data: 3 observations and 250 variables (genes)
set.seed(724)
xN <- matrix(rnorm(750), nrow = 3, ncol = 250)

## Fitting Wilma
fit  <- wilma(leukemia.x, leukemia.y, noc = 3, trace = 1)

## Working with the output
fit
summary(fit)
plot(fit)
fitted(fit)

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
\keyword{cluster}
