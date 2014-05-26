\name{run.BayesianLR.test}
\alias{run.BayesianLR.test}
\title{The actual function that does the running (BayesianLR.test is now a wrapper that runs this with ever-looser tolerances)}
\usage{
  run.BayesianLR.test(truePos, totalDzPos, trueNeg,
    totalDzNeg, R = 5 * 10^4, verbose = FALSE,
    parameters = list(shrink = 5, tol = 5e-04, nEach = 80),
    ...)
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

  \item{\dots}{Arguments to pass along to boot.ci for the
  BCa confidence intervals.}
}
\value{
  An object of class lrtest.
}
\description{
  The actual function that does the running
  (BayesianLR.test is now a wrapper that runs this with
  ever-looser tolerances)
}

