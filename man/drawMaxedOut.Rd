\name{drawMaxedOut}
\alias{drawMaxedOut}
\title{Internal function to draw a set of sensitivities or specificities
This is intended for the case where testPos == totalDzPos or testNeg == totalDzNeg.}
\usage{
  drawMaxedOut(n, R, verbose,
    parameters = list(shrink = 5, tol = 5e-04, nEach = 80))
}
\arguments{
  \item{n}{The total number of positives/negatives in the
  population.}

  \item{R}{is the number of replications in each round of
  the bootstrap (has been tested at 50,000 or greater).}

  \item{verbose}{Whether to display internal operations as
  they happen.}

  \item{parameters}{List of control parameters (shrink,
  tol, nEach) for sequential grid search.}
}
\description{
  Internal function to draw a set of sensitivities or
  specificities This is intended for the case where testPos
  == totalDzPos or testNeg == totalDzNeg.
}

