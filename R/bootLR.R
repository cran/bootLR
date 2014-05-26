# ----- Single-valued calculations ----- #

#' Compute sensitivity, specificity, positive likelihood ratio, negative likelihood ratio for a single 2x2 table
#' @param truePos The number of true positive tests.
#' @param totalDzPos The total number of positives ("sick") in the population.
#' @param trueNeg The number of true negatives in the population.
#' @param totalDzNeg The total number of negatives ("well") in the population.
#' @return A one-row matrix containing sensitivity, specificity, posLR, negLR results.
#' @references Deeks JJ, Altman DG. BMJ. 2004 July 17; 329(7458): 168-169.
#' @examples
#' \dontrun{
#' confusionStatistics( 25, 50, 45, 75 )
#' }
confusionStatistics <- function( truePos, totalDzPos, trueNeg, totalDzNeg ) {
  n <- length(truePos)
  res <- matrix( NA, ncol=4, nrow=n )
  colnames(res) <- c("sens","spec","posLR","negLR")
  res[,"sens"] <- truePos / totalDzPos
  res[,"spec"] <- trueNeg / totalDzNeg
  res[,"posLR"] <- res[,"sens"] / ( 1 - res[,"spec"] )
  res[,"negLR"] <- ( 1 - res[,"sens"] ) / res[,"spec"]
  res
}

# ----- Optimization tools ----- #

#' Find the lowest population probability whose median is consistently one
#' This is the lowest estimate for Sens that is consistently (over 5 runs) most likely to yield a sample estimate that is all 1's (e.g. 100/100, etc.).
#' @param pr Probability input.
#' @param size Number of trials.
#' @param R number of bootstrap replications.
#' @param nConsistentRuns Number of runs that all have to be identical to return TRUE.
#' @param warn Warn if searching outside of the range c(0,1).
#' @return Boolean of length one (TRUE or FALSE).
#' @examples
#' \dontrun{
#' prs <- seq(.990,.995,.0001)
#' bools <- sapply( prs, medianConsistentlyOne, size=truePos, R=R )
#' data.frame( prs, bools )
#' }
medianConsistentlyOne <- function(pr, size, R, nConsistentRuns=5, warn=TRUE) {
  if( 0 > pr | pr > 1) {
    if(warn)  warning("Searching probabilities outside of 0,1. Returning FALSE.")
    return( FALSE )
  } else {
    reps <- replicate( nConsistentRuns, median( rbinom(R, size=size, prob=pr) ) )
    return( all( reps==size ) )
  }
}

#' Optimize a function returning a single numeric value subject to a boolean constraint
#' Utilizes a naive recursive grid search.
#' @param f Function to be minimized: takes a single numeric value and returns a single numeric value.
#' @param constraint Function of a single variable returning a single boolean value (must be TRUE to be at the optimum).
#' @param bounds A numeric vector of length two which are the upper and lower bounds of the input to try.
#' @param nEach Number of points n each round of grid searching to use.
#' @param shrink Factor indicating how much (1/shrink) to narrow the search width by each round; highly recommended that shrink is at least half the size of nEach.
#' @param tol The tolerance (epsilon).
#' @param verbose Whether to display verbose output.
#' @param \dots Arguments to pass along to constraint.
#' @return The optimized input value (numeric).
sequentialGridSearch <- function( f, constraint, bounds, nEach=40, shrink=10, tol=.Machine$double.eps ^ 0.5, verbose=FALSE, ... ) {
  #! The alabama package or similar might be a better way of doing this in the future.
  if(verbose) cat("Grid searching between",bounds[1],"and",bounds[2],"\n")
  x <- seq( from=bounds[1], to=bounds[2], length.out=nEach )
  fx <- f( x )
  if( any(is.na(fx)) ) stop("NAs produced while evaluating f")
  cx <- constraint( x, ... )
  if( any(is.na(cx)) ) stop("NAs produced while evaluating constraint")
  if( !any(cx) ) stop("No value found while searching between ",bounds[1], " and ",bounds[2],". Try setting a looser tolerance, a lower shrinkage value, or a higher number for nEach.\n")
  newVal <- x[ which( fx==min( fx[cx] ) ) ] #! Not very efficient
  if(verbose) cat("Newval found as:", newVal,", producing value",f(newVal),"\n")
  if( exists("lastVal") && abs(f(newVal)-f(lastVal))<tol ) { # Successive rounds within tolerance
    if(verbose) cat("Successive rounds within tolerance.  Success!\n")
    return(newVal)
  } else if( sum(cx) > 1 && abs( f(newVal) - min(fx[cx & x!=newVal]) )<tol ) { # Two values this round are within tolerance
    if(verbose) cat("Two values this round within tolerance.  Success!\n")
    return(newVal)
  } else { # No values within tolerance, but at least one value still works--keep recursing!
    if(verbose) cat("No values within tolerance.  Recursing.\n")
    newHalfRange <- ( abs(diff(range(bounds)))/shrink ) / 2
    newBounds <- c( newVal - newHalfRange, newVal + newHalfRange ) # New bounds with narrower range, centered on crude optimum so far
    lastVal <- newVal
    return( sequentialGridSearch( f=f, constraint=constraint, bounds=newBounds, nEach=nEach, shrink=shrink, tol=tol, verbose=verbose, ... ) ) 
  }
}

# ----- Main function and its helpers ----- #

#' Compute the (positive/negative) likelihood ratio with appropriate, bootstrapped confidence intervals
#' 
#' Compute the (positive/negative) likelihood ratio with appropriate, bootstrapped confidence intervals. 
#' A standard bootstrapping approach is used for sensitivity and specificity, results are combined, and 
#' then 95% CIs are determined. 
#' For the case where sensitivity or specificity equals zero or one, an appropriate bootstrap sample is generated 
#' and then used in subsequent computations.  
#' 
#' If the denominator is 0, calculations are inverted until the final result.
#' 
#' @param truePos The number of true positive tests.
#' @param totalDzPos The total number of positives ("sick") in the population.
#' @param trueNeg The number of true negatives in the population.
#' @param totalDzNeg The total number of negatives ("well") in the population.
#' @param R is the number of replications in each round of the bootstrap (has been tested at 50,000 or greater).
#' @param verbose Whether to display internal operations as they happen.
#' @param parameters List of control parameters (shrink, tol, nEach) for sequential grid search.
#' @param maxTries Each time a run fails, BayesianLR.test will back off on the parameters and try again. maxTries specifies the number of times to try before giving up.  If you can't get it to converge, try setting this higher.
#' @param \dots Arguments to pass along to boot.ci for the BCa confidence intervals.
#' @return An object of class lrtest.
#' @export BayesianLR.test
#' @imports boot
#' @examples
#' blrt <- BayesianLR.test( truePos=100, totalDzPos=100, trueNeg=60, totalDzNeg=100 )
#' blrt
#' summary(blrt)
#' \dontrun{
#' BayesianLR.test( truePos=98, totalDzPos=100, trueNeg=60, totalDzNeg=100 )
#' BayesianLR.test( truePos=60, totalDzPos=100, trueNeg=100, totalDzNeg=100 )
#' BayesianLR.test( truePos=60, totalDzPos=100, trueNeg=99, totalDzNeg=100 )
#' # Note the argument names are not necessary if you specify them in the proper order:
#' BayesianLR.test( 60, 100, 50, 50 ) 
#' # You can specify R= to increase the number of bootstrap replications
#' BayesianLR.test( 60, 100, 50, 50, R=10000 ) 
#' }
#' @note This algorithm utilizes a sequential grid search.  You'll either need a fast computer or substantial patience for certain combinations of inputs.
BayesianLR.test <- function( truePos, totalDzPos, trueNeg, totalDzNeg, R=5*10^4, verbose=FALSE, parameters=list(shrink=5,tol=.0005,nEach=80), maxTries = 20, ... ) {
  convergeFailText <- "try setting a looser tolerance, a lower shrinkage value, or a higher number for neach" # Error text that indicates a failure of convergence
  res <- structure(NULL,class="try-error",condition=convergeFailText)
  tries <- 1
  while( class(res) == "try-error"  &  tries < maxTries ) {
    if( verbose & tries > 1 )  message("Failed to reach convergence in trial number ", tries-1, ".\nRunning trial number ", tries, " to see if we can reach convergence. New parameters: \nShrink ", parameters$shrink, "\nTolerance ", parameters$tol, "\nnEach ", parameters$nEach,"\n" )
    res <- try( run.BayesianLR.test(truePos, totalDzPos, trueNeg, totalDzNeg, R, verbose, parameters) )
    if( class(res) == "try-error"  &&  !grepl( convergeFailText, tolower( as.character( attributes(res)$condition ) ) ) )  stop( as.character( attributes(res)$condition ) )
    parameters$tol <- ifelse( parameters$tol > .001, parameters$tol, .001 )
    parameters$shrink <- (parameters$shrink - 1) * .65 + 1
    parameters$nEach <- floor( parameters$nEach * 1.3 )
    tries <- tries + 1    
  }
  res
}

#' The actual function that does the running (BayesianLR.test is now a wrapper that runs this with ever-looser tolerances)
#' @param truePos The number of true positive tests.
#' @param totalDzPos The total number of positives ("sick") in the population.
#' @param trueNeg The number of true negatives in the population.
#' @param totalDzNeg The total number of negatives ("well") in the population.
#' @param R is the number of replications in each round of the bootstrap (has been tested at 50,000 or greater).
#' @param verbose Whether to display internal operations as they happen.
#' @param parameters List of control parameters (shrink, tol, nEach) for sequential grid search.
#' @param \dots Arguments to pass along to boot.ci for the BCa confidence intervals.
#' @return An object of class lrtest.
run.BayesianLR.test <- function( truePos, totalDzPos, trueNeg, totalDzNeg, R=5*10^4, verbose=FALSE, parameters=list(shrink=5,tol=.0005,nEach=80), ... ) {
  # -- Check inputs -- #
  if( R < 5*10^4 ) warning("Setting the number of bootstrap replications to a number lower than 50,000 may lead to unstable results")
  if( totalDzPos == 0 | totalDzNeg == 0 ) stop("This package may seem like magic, but not even magic will solve your problem (totalDzPos or totalDzNeg = 0).")
  if( trueNeg > totalDzNeg | truePos > totalDzPos ) stop("You cannot have more test positive/negative than you have total positive/negative.")
  
  # -- Bootstrap sensitivity and specificity -- #
  cs <- confusionStatistics( truePos=truePos, totalDzPos=totalDzPos, trueNeg=trueNeg, totalDzNeg=totalDzNeg )
  csExact <- cs # store the actual confusion statistics, since we will use the lprb as a proxy for them at various points but we still want to report the real numbers at the end
  
  bootmean <- function(x,i)  mean(x[i])
  
  if( truePos == totalDzPos ) {
    sensb <- drawMaxedOut( n=totalDzPos, R=R, verbose=verbose )
    cs[,"sens"] <- attr(sensb,"lprb")
  } else {
    sensb <- boot::boot(
      rep( 1:0, c( truePos, totalDzPos-truePos ) ), 
      bootmean, 
      R=R
    )$t
  }
  
  if( trueNeg == totalDzNeg ) {
    specb <- drawMaxedOut( n=totalDzNeg, R=R, verbose=verbose, parameters=parameters )
    cs[,"spec"] <- attr(specb,"lprb")
  } else {
    specb <- boot::boot(
      rep( 1:0, c( trueNeg, totalDzNeg-trueNeg ) ), 
      bootmean, 
      R=R
    )$t
  }
  
  # -- Compute pos/neg LRs and their BCa confidence intervals -- #
  negLR <- unname( ( 1 - cs[,"sens"] ) / cs[,"spec"]  )
  negLRexact <- unname( ( 1-csExact[,"sens"] ) / csExact[,"spec"] ) # Could also just use csExact[,"negLR"] here, but not in the line above it
  if( all( specb != 0L ) ) {
    negLR.ci <- bca( ( 1 - sensb) / specb, negLR, ... )$bca[4:5]
  } else {
    negLR.ci <- 1/bca( specb / ( 1 - sensb), 1/negLR, ... )$bca[c(4,5)]
  }
  posLR <- unname( cs[,"sens"] / ( 1 - cs[,"spec"] ) )
  posLRexact <- unname( csExact[,"sens"] / ( 1 - csExact[,"spec"] ) )
  if( all( specb != 1L ) ) {
    posLR.ci <- bca( sensb / ( 1 - specb ), posLR, ... )$bca[4:5]
  } else {
    posLR.ci <- 1/bca( ( 1 - specb ) / sensb, 1/posLR, ... )$bca[c(5,4)] # Reversed because the order inverts when you take the reciprocal
  }
  
  # -- Return lrtest object -- #
  structure( list(
    negLR = negLRexact,
    negLR.ci = negLR.ci,
    posLR = posLRexact,
    posLR.ci = posLR.ci,
    inputs = structure( c( truePos, totalDzPos, trueNeg, totalDzNeg ), names=c("truePos","totalDzPos","trueNeg","totalDzNeg") ),
    statistics = cs[ , c("sens","spec") ]
  ), 
  class = "lrtest",
  ci.type = "BCa",
  ci.width = .95
  )
}

#' Internal function to draw a set of sensitivities or specificities
#' This is intended for the case where testPos == totalDzPos or testNeg == totalDzNeg.
#' @param n The total number of positives/negatives in the population.
#' @param R is the number of replications in each round of the bootstrap (has been tested at 50,000 or greater).
#' @param verbose Whether to display internal operations as they happen.
#' @param parameters List of control parameters (shrink, tol, nEach) for sequential grid search.
drawMaxedOut <- function( n, R, verbose, parameters=list(shrink=5,tol=.0005,nEach=80) ) {
  lprb <- sequentialGridSearch( # lowest probability that consistently produces 1's 
    f=identity, # We just want to minimize pr
    constraint=function(probs,...) vapply( probs, FUN=medianConsistentlyOne, FUN.VALUE=NA, ... ),
    bounds=c(0,1), 
    verbose=verbose,
    size=n, R=R, warn=FALSE,
    shrink=parameters$shrink,
    tol=parameters$tol,
    nEach=parameters$nEach
  )
  res <- rbinom(R, size=n, prob=lprb)/n
  attr( res, "lprb" ) <- lprb
  res
}


#' Internal function to analyze LR bootstrap finding median, and standard and
#' BCa percentile 95% CIs.
#' To obtain bca CI on a non-boot result, use a dummy boot.
#' and replace t and t0 with the results of interest.
#' @param t The vector to obtain a BCa bootstrap for (e.g. nlr).
#' @param t0 The central value of the vector (e.g. the ).
#' @param \dots Pass-alongs to boot.ci.
bca <- function( t, t0, ... ) {
  R <- length(t)
  dummy <- rep(1:0,c(5,5)) # Doesn't matter what values are given here, since we're replacing them
  dummyb <- boot::boot(dummy, function(x,i) 1, R=R)
  dummyb$t <- matrix(t,ncol=1)
  dummyb$t0 <- t0
  boot::boot.ci(dummyb, t0=dummyb$t0, t=dummyb$t, type=c("perc", "bca"), ...)
}

# ----- Functions to display the resulting lrtest object ----- #

#' Prints results from the BayesianLR.test
#' As is typical for R, this is run automatically when you type in an object name, and is typically not run directly by the end-user.
#' @param x The lrtest object created by BayesianLR.test.
#' @param \dots Pass-alongs (currently ignored).
#' @return Returns x unaltered.
#' @method print lrtest
#' @S3method print lrtest
#' @export print.lrtest
print.lrtest <- function( x, ... ) {
  digits <- 3 # Number of digits to round to for display purposes
  cat("\n")
  cat("Likelihood ratio test of a 2x2 table")
  cat("\n\n")
  cat("data:\n")
  print(x$inputs)
  cat( paste0( "Positive LR: ", round(x$posLR,digits), " (", round(x$posLR.ci[1],digits), " - ", round(x$posLR.ci[2],digits), ")\n" ) )
  cat( paste0( "Negative LR: ", round(x$negLR,digits), " (", round(x$negLR.ci[1],digits), " - ", round(x$negLR.ci[2],digits), ")\n" ) )
  cat( paste0( attr(x,"ci.width")*100, "% confidence intervals computed via ", attr(x,"ci.type"), " bootstrapping.\n" ) )
  cat( "Note: This procedure depends on repeated random sampling.  As such it is subject to some variability in results.\n  Variability is minimized by large numbers of replications (generally 50,000) [and averaging 5 repeated results],\n but with small sample sizes or sensitivity or specificity near 0 or 1, variability becomes more pronounced.\n  This is not an error, it is a function of the nature of the procedure." )
  invisible(x)
}