\name{BayesianLR.test}
\alias{BayesianLR.test}
\title{Compute the (positive/negative) likelihood ratio with appropriate, bootstrapped confidence intervals}
\usage{
  BayesianLR.test(truePos, totalDzPos, trueNeg, totalDzNeg,
    R = 5 * 10^4, verbose = FALSE,
    parameters = list(shrink = 5, tol = 5e-04, nEach = 80),
    maxTries = 20, ...)
}
\arguments{
  \item{truePos}{The number of true positive tests.}

  \item{totalDzPos}{The total number of positives ("sick")
  in the population.}

  \item{trueNeg}{The number of true negatives in the
  population.}

  \item{totalDzNeg}{The total number of negatives ("well")
  in the population.}

  \item{R}{is the number of replications in each round of
  the bootstrap (has been tested at 50,000 or greater).}

  \item{verbose}{Whether to display internal operations as
  they happen.}

  \item{parameters}{List of control parameters (shrink,
  tol, nEach) for sequential grid search.}

  \item{maxTries}{Each time a run fails, BayesianLR.test
  will back off on the parameters and try again. maxTries
  specifies the number of times to try before giving up.
  If you can't get it to converge, try setting this
  higher.}

  \item{\dots}{Arguments to pass along to boot.ci for the
  BCa confidence intervals.}
}
\value{
  An object of class lrtest.
}
\description{
  Compute the (positive/negative) likelihood ratio with
  appropriate, bootstrapped confidence intervals. A
  standard bootstrapping approach is used for sensitivity
  and specificity, results are combined, and then 95% CIs
  are determined. For the case where sensitivity or
  specificity equals zero or one, an appropriate bootstrap
  sample is generated and then used in subsequent
  computations.
}
\details{
  If the denominator is 0, calculations are inverted until
  the final result.
}
\note{
  This algorithm utilizes a sequential grid search.  You'll
  either need a fast computer or substantial patience for
  certain combinations of inputs.
}
\examples{
blrt <- BayesianLR.test( truePos=100, totalDzPos=100, trueNeg=60, totalDzNeg=100 )
blrt
summary(blrt)
\dontrun{
BayesianLR.test( truePos=98, totalDzPos=100, trueNeg=60, totalDzNeg=100 )
BayesianLR.test( truePos=60, totalDzPos=100, trueNeg=100, totalDzNeg=100 )
BayesianLR.test( truePos=60, totalDzPos=100, trueNeg=99, totalDzNeg=100 )
# Note the argument names are not necessary if you specify them in the proper order:
BayesianLR.test( 60, 100, 50, 50 )
# You can specify R= to increase the number of bootstrap replications
BayesianLR.test( 60, 100, 50, 50, R=10000 )
}
}

