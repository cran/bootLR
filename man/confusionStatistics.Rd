\name{confusionStatistics}
\alias{confusionStatistics}
\title{Compute sensitivity, specificity, positive likelihood ratio, negative likelihood ratio for a single 2x2 table}
\usage{
  confusionStatistics(truePos, totalDzPos, trueNeg,
    totalDzNeg)
}
\arguments{
  \item{truePos}{The number of true positive tests.}

  \item{totalDzPos}{The total number of positives ("sick")
  in the population.}

  \item{trueNeg}{The number of true negatives in the
  population.}

  \item{totalDzNeg}{The total number of negatives ("well")
  in the population.}
}
\value{
  A one-row matrix containing sensitivity, specificity,
  posLR, negLR results.
}
\description{
  Compute sensitivity, specificity, positive likelihood
  ratio, negative likelihood ratio for a single 2x2 table
}
\examples{
\dontrun{
confusionStatistics( 25, 50, 45, 75 )
}
}
\references{
  Deeks JJ, Altman DG. BMJ. 2004 July 17; 329(7458):
  168-169.
}

