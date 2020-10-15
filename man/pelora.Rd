\name{pelora}
\alias{pelora}
\title{Supervised Grouping of Predictor Variables}
\description{Performs selection and supervised grouping of predictor
  variables in large (microarray gene expression) datasets, with an
  option for simultaneous classification. Works in a greedy forward
  strategy and optimizes the binomial log-likelihood, based on estimated
  conditional probabilities from penalized logistic regression analysis.}
\usage{
pelora(x, y, u = NULL, noc = 10, lambda = 1/32, flip = "pm",
       standardize = TRUE, trace = 1)
}

\arguments{
  \item{x}{Numeric matrix of explanatory variables (\eqn{p} variables in
    columns, \eqn{n} cases in rows). For example, these can be
    microarray gene expression data which should be grouped.}
  \item{y}{Numeric vector of length \eqn{n} containing the class labels
    of the individuals. These labels have to be coded by 0 and 1.}
  \item{u}{Numeric matrix of additional (clinical) explanatory variables
    (\eqn{m} variables in columns, \eqn{n} cases in rows) that are used
    in the (penalized logistic regression) prediction model, but neither
    grouped nor averaged. For example, these can be 'traditional'
    clinical variables.}
  \item{noc}{Integer, the number of clusters that should be searched for
    on the data.}
  \item{lambda}{Real, defaults to 1/32. Rescaled penalty parameter that
    should be in \eqn{[0,1]}.}
  \item{flip}{Character string, describing a method how the \code{x}
    (gene expression) matrix should be sign-flipped. Possible are
    \code{"pm"} (the default) where the sign for each variable is
    determined upon its entering into the group, \code{"cor"} where the
    sign for each variable is determined a priori as the sign of the
    empirical correlation of that variable with the \code{y}-vector, and
    \code{"none"} where no sign-flipping is carried out.}
  \item{standardize}{Logical, defaults to \code{TRUE}. Is indicating
    whether the predictor variables (genes) should be standardized to
    zero mean and unit variance.}
  \item{trace}{Integer >= 0; when positive, the output of the internal
    loops is provided; \code{trace >= 2} provides output even from the
    internal C routines.}
}
\value{
  \code{pelora} returns an object of class "pelora". The functions
  \code{print} and \code{summary} are used to obtain an overview of the
  variables (genes) that have been selected and the groups that have
  been formed. The function \code{plot} yields a two-dimensional
  projection into the space of the first two group centroids that
  \code{pelora} found. The generic function \code{fitted} returns
  the fitted values, these are the cluster representatives. \code{coef}
  returns the penalized logistic regression coefficients \eqn{\theta_j}
  for each of the predictors. Finally, \code{predict} is used for
  classifying test data with Pelora's internal penalized logistic
  regression classifier on the basis of the (gene) groups that have been
  found.

  An object of class "pelora" is a list containing:
  \item{genes}{A list of length \code{noc}, containing integer vectors
    consisting of the indices (column numbers) of the variables (genes)
    that have been clustered.}
  \item{values}{A numerical matrix with dimension \eqn{n \times
      \code{noc}}, containing the fitted values, i.e. the group
    centroids \eqn{\tilde{x}_j}.}
  \item{y}{Numeric vector of length \eqn{n} containing the class labels
    of the individuals. These labels are coded by 0 and 1.}
  \item{steps}{Numerical vector of length \code{noc}, showing the number
    of forward/backward cycles in the fitting process of each cluster.}
  \item{lambda}{The rescaled penalty parameter.}
  \item{noc}{The number of clusters that has been searched for on the data.}
  \item{px}{The number of columns (genes) in the \code{x}-matrix.}
  \item{flip}{The method that has been chosen for sign-flipping the
    \code{x}-matrix.}
  \item{var.type}{A factor with \code{noc} entries, describing whether
    the \eqn{j}th predictor is a group of predictors (genes) or a single
    (clinical) predictor variable.}
  \item{crit}{A list of length \code{noc}, containing numerical vectors
    that provide information about the development of the grouping
    criterion during the clustering.}
  \item{signs}{Numerical vector of length \eqn{p}, saying whether the
    \eqn{i}th variable (gene) should be sign-flipped (-1) or not (+1).}
  \item{samp.names}{The names of the samples (rows) in the
    \code{x}-matrix.}
  \item{gene.names}{The names of the variables (columns) in the
    \code{x}-matrix.}
  \item{call}{The function call.}
}
\references{
  Marcel Dettling (2003)
  \emph{Finding Predictive Gene Groups from Microarray Data}, see
  \url{https://stat.ethz.ch/~dettling/supervised.html}

  Marcel Dettling and Peter Bühlmann (2002).
  Supervised Clustering of Genes.
  \emph{Genome Biology}, \bold{3}(12): research0069.1-0069.15,
  \doi{10.1186/gb-2002-3-12-research0069}.

  Marcel Dettling and Peter Bühlmann (2004).
  Finding Predictive Gene Groups from Microarray Data.
  \emph{Journal of Multivariate Analysis} \bold{90}, 106--131,
  \doi{10.1016/j.jmva.2004.02.012}
}
\author{Marcel Dettling, \email{dettling@stat.math.ethz.ch}}
\seealso{\code{\link{wilma}} for another supervised clustering technique.}
\examples{
## Working with a "real" microarray dataset
data(leukemia, package="supclust")

## Generating random test data: 3 observations and 250 variables (genes)
set.seed(724)
xN <- matrix(rnorm(750), nrow = 3, ncol = 250)

## Fitting Pelora
fit <- pelora(leukemia.x, leukemia.y, noc = 3)

## Working with the output
fit
summary(fit)
plot(fit)
fitted(fit)
coef(fit)

## Fitted values and class probabilities for the training data
predict(fit, type = "cla")
predict(fit, type = "prob")

## Predicting fitted values and class labels for the random test data
predict(fit, newdata = xN)
predict(fit, newdata = xN, type = "cla", noc = c(1,2,3))
predict(fit, newdata = xN, type = "pro", noc = c(1,3))

## Fitting Pelora such that the first 70 variables (genes) are not grouped
fit <- pelora(leukemia.x[, -(1:70)], leukemia.y, leukemia.x[,1:70])

## Working with the output
fit
summary(fit)
plot(fit)
fitted(fit)
coef(fit)

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


