# ----- Single-valued calculations ----- #

#' Compute sensitivity, specificity, positive likelihood ratio, negative likelihood ratio for a single 2x2 table
#' @param truePos The number of true positive tests.
#' @param totalDzPos The total number of positives ("sick") in the population.
#' @param trueNeg The number of true negatives in the population.
#' @param totalDzNeg The total number of negatives ("well") in the population.
#' @return A one-row matrix containing sensitivity, specificity, posLR, negLR results.
#' @references Deeks JJ, Altman DG. BMJ. 2004 July 17; 329(7458): 168-169.
#' @export confusionStatistics
#' @examples
#' \dontrun{
#' confusionStatistics( 25, 50, 45, 75 )
#' }
confusionStatistics <- function( truePos, totalDzPos, trueNeg, totalDzNeg ) {
  n <- length(truePos)
  if( length(truePos) != length(totalDzPos) | length(truePos) != length(trueNeg) | length(truePos) != length(totalDzNeg) )  stop("Lengths of the input vectors must match.")
  res <- matrix( NA, ncol=4, nrow=n )
  colnames(res) <- c("sens","spec","posLR","negLR")
  res[,"sens"] <- truePos / totalDzPos
  res[,"spec"] <- trueNeg / totalDzNeg
  res[,"posLR"] <- res[,"sens"] / ( 1 - res[,"spec"] )
  res[,"negLR"] <- ( 1 - res[,"sens"] ) / res[,"spec"]
  res
}

#' Compute values and confidence intervals for sensitivity, specificity, positive likelihood ratio, negative likelihood ratio for a single 2x2 table
#' @param truePos The number of true positive tests.
#' @param totalDzPos The total number of positives ("sick") in the population.
#' @param trueNeg The number of true negatives in the population.
#' @param totalDzNeg The total number of negatives ("well") in the population.
#' @param calcLRCI Method to use to calculate the LR CI: "BayesianLR.test" "none" or "analytic"
#' @param \dots Arguments to pass to Bayesian.LRtest. 
#' @param alpha The alpha for the width of the confidence interval (defaults to alpha = 0.05 for a 95 percent CI)
#' @param binomMethod The method to be passed to binom.confint to calculate confidence intervals of proportions (sensitivity, etc.).  See help("binom.confint") and the Newcombe article referenced below.
#' @return A matrix containing sensitivity, specificity, posLR, negLR results and their confidence intervals
#' @references Deeks JJ, Altman DG. BMJ. 2004 July 17; 329(7458): 168-169. Newcombe RG. Statist Med. 1998; 17(857-872).
#' @export diagCI
#' @examples
#' \dontrun{
#' diagCI( 25, 50, 45, 75 )
#' diagCI( truePos = c(25, 30), totalDzPos = c( 50, 55 ), trueNeg = c(5, 35), totalDzNeg = c(60,65) )
#' }
diagCI <- function( truePos, totalDzPos, trueNeg, totalDzNeg, calcLRCI = "BayesianLR.test", alpha = 0.05, binomMethod = "wilson", ... ) {
  # Validate inputs
  if( !calcLRCI %in% c("BayesianLR.test", "none", "analytic") )  stop( "calcLRCI is not a valid option" )
  if( alpha <= 0 | alpha >= 1 )  stop( "alpha must be between 0 and 1" )
  if( binomMethod == "all" | length( binomMethod ) != 1 )  stop( "Only one method allowed" )
  if( length(truePos) != length(totalDzPos) | length(truePos) != length(trueNeg) | length(truePos) != length(totalDzNeg) )  stop("Lengths of the input vectors must match.")
  # Other values needed
  conf.level <- 1 - alpha
  Z <- stats::qnorm( 1 - alpha / 2 )
  n <- length(truePos)
  # Calculate the rest of the 2x2 table
  falsePos <- totalDzNeg - trueNeg # "b"
  falseNeg <- totalDzPos - truePos # "c"
  # Define a function to return the proportion and its confidence interval
  calcPropCI <- function( x, n, binomMethod, conf.level, propname = "prop" ) {
    dat <- data.frame( prop = x / n )
    binomCI <- binom::binom.confint( x, n, methods = binomMethod, conf.level = conf.level )
    dat$prop_LB <- binomCI$lower
    dat$prop_UB <- binomCI$upper
    colnames( dat ) <- sub( "^prop", propname, colnames( dat ) )
    dat
  }
  # Run that function across all statistics
  statInputs <- list( 
    list( propname = "sens", x = truePos, n = totalDzPos ),
    list( propname = "spec", x = trueNeg, n = totalDzNeg ),
    list( propname = "PPV", x = truePos, n = truePos + falsePos ),
    list( propname = "NPV", x = trueNeg, n = trueNeg + falseNeg )
  )
  resList <- lapply( statInputs, FUN = function( l, binomMethod, conf.level ) {
    calcPropCI( l$x, l$n, binomMethod = binomMethod, conf.level = conf.level, propname = l$propname )
  }, binomMethod = binomMethod, conf.level = conf.level )
  res <- do.call( cbind, resList )
  # Positive/Negative LR point estimates
  cs <- mapply(
    FUN = confusionStatistics,
    truePos = truePos, totalDzPos = totalDzPos, trueNeg = trueNeg, totalDzNeg = totalDzNeg,
    SIMPLIFY = FALSE
  )
  csDF <- t( vapply( cs, FUN = function(x) unlist(x[ 1, grep( "LR", colnames(x) ) ]), FUN.VALUE = rep(NA_real_,2) ) ) 
  na <- rep( NA_real_, length( truePos ) )
  res <- cbind( res, data.frame( posLR = na, posLR_LB = na, posLR_UB = na, negLR = na, negLR_LB = na, negLR_UB = na ) )
  res[,"posLR"] <- csDF[,"posLR"]
  res[,"negLR"] <- csDF[,"negLR"]
  # Positive/Negative LR confidence intervals
  if( calcLRCI == "BayesianLR.test" ) {
    blrt <- mapply(
      FUN = BayesianLR.test,
      truePos = truePos, totalDzPos = totalDzPos, trueNeg = trueNeg, totalDzNeg = totalDzNeg,
      SIMPLIFY = FALSE,
      ...
    )
    blrtDF <- t( vapply( blrt, FUN = function(x) unlist(x[ grep( "LR", names(x) ) ]), FUN.VALUE = rep(NA_real_,6) ) ) 
    res[,"posLR"] <- blrtDF[,"posLR"]
    res[,"posLR_LB"] <- apply( blrtDF[ , grep( "posLR.ci", colnames( blrtDF ) ), drop = FALSE ], 1, min ) 
    res[,"posLR_UB"] <- apply( blrtDF[ , grep( "posLR.ci", colnames( blrtDF ) ), drop = FALSE ], 1, max ) 
    res[,"negLR"] <- blrtDF[,"negLR"]
    res[,"negLR_LB"] <- apply( blrtDF[ , grep( "negLR.ci", colnames( blrtDF ) ), drop = FALSE ], 1, min ) 
    res[,"negLR_UB"] <- apply( blrtDF[ , grep( "negLR.ci", colnames( blrtDF ) ), drop = FALSE ], 1, max ) 
  } else if( calcLRCI == "analytic" ) {
    # Analytic solution to LR CI
    posLRSE <- sqrt( (1/truePos) - (1/(truePos+falseNeg)) + (1/falsePos) - (1/(falsePos+trueNeg))  )
    res[,"posLR_LB"] <- exp( log( res[,"posLR"] ) - Z * posLRSE )
    res[,"posLR_UB"] <- exp( log( res[,"posLR"] ) + Z * posLRSE )
    negLRSE <- sqrt( (1/falseNeg) - (1/(truePos+falseNeg)) + (1/trueNeg) - (1/(falsePos+trueNeg)) )
    res[,"negLR_LB"] <- exp( log( res[,"negLR"] ) - Z * negLRSE )
    res[,"negLR_UB"] <- exp( log( res[,"negLR"] ) + Z * negLRSE )
  }
  # Add the input data
  res <- cbind( truePos, totalDzPos, trueNeg, totalDzNeg, res )
  # Return
  structure( res, class = c("diagCI","data.frame"), calcLRCI = calcLRCI, alpha = alpha, conf.level = conf.level, binomMethod = binomMethod )
}

#' Prints results from diagCI
#' As is typical for R, this is run automatically when you type in an object name, and is typically not run directly by the end-user.
#' @param x The diagCI object created by diagCI()
#' @param digits Number of digits to round to
#' @param \dots Pass-alongs (currently ignored).
#' @return Returns x unaltered.
#' @method print diagCI 
#' @export
#' @examples
#' \dontrun{
#' diagCI( 25, 50, 45, 75 )
#' }
print.diagCI <- function( x, digits = 3, ... ) {
  apply( x, 1, function( rw ) {
    # Print the 2x2 table
    twoByTwo <- matrix( data = c( rw["truePos"], rw["totalDzPos"] - rw["truePos"], rw["totalDzNeg"] - rw["trueNeg"], rw["trueNeg"] ), nrow = 2, ncol = 2, dimnames = list( c("Test Positive", "Test Negative"), c("Condition Positive", "Condition Negative") ) )
    print( twoByTwo )
    # Print the diagnostic statistics
    rw <- rw[ seq( -1, -4, -1 ) ]
    rwMat <- matrix( rw, nrow = 3 )
    rwMat <- round( rwMat, digits )
    rwMat <- rbind( names( rw )[ seq( 1, length( rw ), 3 ) ], rwMat )
    statistics <- apply( rwMat, 2, function( d ) {
      paste0( d[1], " ", d[2], " ( ", d[3], " - ", d[4], " )" )
    } )
    for (s in statistics) { print( s ) }
  } )
  invisible( x )
}

# ----- Optimization tools ----- #

#' Find the lowest population probability whose median is consistently one
#' This is the lowest estimate for Sens that is consistently (over 5 runs) most likely to yield a sample estimate that is all 1's (e.g. 100/100, etc.).
#' @param pr Probability input.
#' @param size Number of trials.
#' @param R number of bootstrap replications.
#' @param nConsistentRuns Number of runs that all have to be identical to return TRUE.
#' @param warn Warn if searching outside of the range c(0,1).
#' @param consistentQuantile Defaults to 0.5 (the median). Change if we want to use a different criterion for consistency than the median
#' @return Boolean of length one (TRUE or FALSE).
#' @examples
#' \dontrun{
#' prs <- seq(.990,.995,.0001)
#' bools <- sapply( prs, medianConsistentlyOne, size=truePos, R=R )
#' data.frame( prs, bools )
#' }
medianConsistentlyOne <- function(pr, size, R, nConsistentRuns=5, warn=TRUE, consistentQuantile = 0.5) {
  if( 0 > pr | pr > 1) {
    if(warn)  warning("Searching probabilities outside of 0,1. Returning FALSE.")
    return( FALSE )
  } else {
    reps <- replicate( nConsistentRuns, stats::quantile( stats::rbinom(R, size=size, prob=pr), probs = consistentQuantile, names = FALSE, type = 8 ) )
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
#' @param R The number of replications in each round of the bootstrap (has been tested at 10,000 or greater).
#' @param nBSave The number of times to re-bootstrap the statistic and then average at the end to obtain confidence intervals (has been tested at 50).
#' @param verbose Whether to display internal operations as they happen.
#' @param parameters List of control parameters (shrink, tol, nEach) for sequential search.
#' @param maxTries Each time a run fails, BayesianLR.test will back off on the parameters and try again. maxTries specifies the number of times to try before giving up.  If you can't get it to converge, try setting this higher.
#' @param ci.width Changing this parameter results in different properties than have been tested and is not recommended. The width of the confidence interval used by boot.ci (not necessarily the same as the width of the CI produced by the algorithm overall)
#' @param consistentQuantile Changing this parameter results in different properties than have been tested and is not recommended. Defaults to 0.5, i.e. the median. Finds the lowest probability for which the random draws are likely to be consistently one, where consistently is defined by this value (i.e. at .5, a simple majority of the time is enough for consistency).
#' @param \dots Arguments to pass along to boot.ci for the BCa confidence intervals.
#' @return An object of class lrtest.
#' @export BayesianLR.test
#' @examples
#' \dontrun{
#' blrt <- BayesianLR.test( truePos=100, totalDzPos=100, trueNeg=60, totalDzNeg=100 )
#' blrt
#' summary(blrt)
#' 
#' BayesianLR.test( truePos=98, totalDzPos=100, trueNeg=60, totalDzNeg=100 )
#' BayesianLR.test( truePos=60, totalDzPos=100, trueNeg=100, totalDzNeg=100 )
#' BayesianLR.test( truePos=60, totalDzPos=100, trueNeg=99, totalDzNeg=100 )
#' 
#' # Note the argument names are not necessary if you specify them in the proper order:
#' BayesianLR.test( 60, 100, 50, 50 ) 
#' 
#' # You can specify R= to increase/decrease the number of bootstrap replications
#' BayesianLR.test( 60, 100, 50, 50, R=10000 ) 
#' 
#' # You can change the number of digits that are printed
#' print.lrtest( BayesianLR.test( 500, 500, 300, 500 ), digits = 4 )
#' 
#' # Or extract the results yourself
#' model.blrt1 <- BayesianLR.test( 500, 500, 300, 500 )
#' unclass( model.blrt1 )
#' model.blrt1$statistics
#' model.blrt1$posLR.ci
#' 
#' # If the model doesn't converge, you can alter the search parameters
#' BayesianLR.test( 500, 500, 300, 500, parameters=list(shrink=4,tol=.001,nEach=150), maxTries = 50 )
#' 
#' ### Statistician-only options
#' # These change the way the model works. 
#' # It is not recommended to alter these, as this will alter the statistical properties of the test
#' # in ways that have not been validated.
#' # Change number of bootstrap replications
#' BayesianLR.test( 500, 500, 300, 500, R = 5*10^4 )
#' # Change number of times to average the confidence interval limits at the end
#' BayesianLR.test( 500, 500, 300, 500, nBSave = 100 )
#' # Change the criteria from median being consistent 0 or 1 to some other quantile
#' BayesianLR.test( 500, 500, 300, 500, consistentQuantile = .53 )
#' }
#' @note You'll either need a fast computer or substantial patience for certain combinations of inputs.
BayesianLR.test <- function( truePos, totalDzPos, trueNeg, totalDzNeg, R=10^4, nBSave=50, verbose=FALSE, parameters=list(shrink=5,tol=.0005,nEach=80), maxTries = 20, ci.width = 0.95, consistentQuantile = 0.5, ... ) {
  # Run the algorithm once, expand criteria if it fails
  convergeFailText <- "try setting a looser tolerance, a lower shrinkage value, or a higher number for neach" # Error text that indicates a failure of convergence
  res <- structure( list() ,class="try-error",condition=convergeFailText)
  tries <- 1
  while( class(res) == "try-error"  &  tries < maxTries ) {
    if( verbose & tries > 1 )  message("Failed to reach convergence in trial number ", tries-1, ".\nRunning trial number ", tries, " to see if we can reach convergence. New parameters: \nShrink ", parameters$shrink, "\nTolerance ", parameters$tol, "\nnEach ", parameters$nEach,"\n" )
    res <- try( run.BayesianLR.test( truePos = truePos, totalDzPos = totalDzPos, trueNeg = trueNeg, totalDzNeg = totalDzNeg, R = R, verbose = verbose, parameters = parameters, ci.width = ci.width, consistentQuantile = consistentQuantile ) )
    if( class(res) == "try-error"  &&  !grepl( convergeFailText, tolower( as.character( attributes(res)$condition ) ) ) )  stop( as.character( attributes(res)$condition ) ) # Unless the error message captured indicates to try again with looser tolerances, fail and return the error message captured from run.BayesianLR.test
    parameters$tol <- ifelse( parameters$tol > .001, parameters$tol, .001 )
    parameters$shrink <- (parameters$shrink - 1) * .65 + 1
    parameters$nEach <- floor( parameters$nEach * 1.3 )
    tries <- tries + 1    
  }
  # Now adjust the posLR and negLR CI's to be the average of nBSave runs of the bootstrap algorithm
  if( verbose )  message("Adjusting CIs by averaging nBSave times.\n")
  bsReps <- replicate( 
    nBSave, 
    unclass( run.BayesianLR.test( truePos = truePos, totalDzPos = totalDzPos, trueNeg = trueNeg, totalDzNeg = totalDzNeg, R = R, verbose = verbose, parameters = parameters, ci.width = ci.width, consistentQuantile = consistentQuantile ) ), 
    simplify=FALSE 
  )
  for( v in c("negLR.ci","posLR.ci") ) {
    res[[ v ]] <- apply( vapply( bsReps, function(x) x[[ v ]], FUN.VALUE = numeric(2) ), 1, mean )
  }
  # Return the result
  res
}

#' The actual function that runs the test. BayesianLR.test is a wrapper that runs this with ever-looser tolerances.
#' @param truePos The number of true positive tests.
#' @param totalDzPos The total number of positives ("sick") in the population.
#' @param trueNeg The number of true negatives in the population.
#' @param totalDzNeg The total number of negatives ("well") in the population.
#' @param R is the number of replications in each round of the bootstrap (has been tested at 10,000 or greater).
#' @param verbose Whether to display internal operations as they happen.
#' @param parameters List of control parameters (shrink, tol, nEach) for sequential grid search.
#' @param ci.width The width of the confidence interval used by boot.ci (not necessarily the same as the width of the CI produced by the algorithm overall; if this parameter is changed, results are not tested)
#' @param consistentQuantile Defaults to 0.5, i.e. the median. Finds the lowest probability for which the random draws are likely to be consistently one, where consistently is defined by this value (i.e. at .5, a simple majority of the time is enough for consistency). Changing this parameter results in different properties than have been tested and is not recommended.
#' @param \dots Arguments to pass along to boot.ci for the BCa confidence intervals.
#' @return An object of class lrtest.
run.BayesianLR.test <- function( truePos, totalDzPos, trueNeg, totalDzNeg, R=10^4, verbose=FALSE, parameters=list(shrink=5,tol=.0005,nEach=80), ci.width = 0.95, consistentQuantile = 0.5, ... ) {
  # -- Check inputs -- #
  if( R < 10^4 ) warning("Setting the number of bootstrap replications to a number lower than 10,000 may lead to unstable results")
  if( totalDzPos == 0 | totalDzNeg == 0 ) stop("This package may seem like magic, but not even magic will solve your problem (totalDzPos or totalDzNeg = 0).")
  if( trueNeg > totalDzNeg | truePos > totalDzPos ) stop("You cannot have more test positive/negative than you have total positive/negative.")
  
  # -- Bootstrap sensitivity and specificity -- #
  cs <- confusionStatistics( truePos=truePos, totalDzPos=totalDzPos, trueNeg=trueNeg, totalDzNeg=totalDzNeg )
  csExact <- cs # store the actual confusion statistics, since we will use the lprb as a proxy for them at various points but we still want to report the real numbers at the end
  
  bootmean <- function(x,i)  mean( x[i] )
  
  if( truePos == totalDzPos ) {
    sensb <- drawMaxedOut( n=totalDzPos, R=R, verbose=verbose, consistentQuantile = consistentQuantile )
    cs[,"sens"] <- attr(sensb,"lprb")
  } else {
    sensb <- boot::boot(
      rep( 1:0, c( truePos, totalDzPos-truePos ) ), 
      bootmean, 
      R=R
    )$t
  }
  
  if( trueNeg == totalDzNeg ) {
    specb <- drawMaxedOut( n=totalDzNeg, R=R, verbose=verbose, parameters=parameters, consistentQuantile )
    cs[,"spec"] <- attr(specb,"lprb")
  } else {
    specb <- boot::boot(
      rep( 1:0, c( trueNeg, totalDzNeg-trueNeg ) ), 
      bootmean, 
      R=R
    )$t
  }
  
  # Test that we don't have any cases where sensb and specb happen to both be 0/1 simultaneously
  if( any(apply( data.frame( sensb, specb ), 1, function(x)  x[1] == 1 & x[2] == 0 )) )
    stop("In some of your draws, sensitivity was 1 at the same time that specificity was 0. This most likely resulted because you had a sensitivity of 1 and a specificity of close to 0, or vice-versa. The algorithm is not designed to handle this case, and may not ever be. Please do not just re-run the algorithm until you no longer receive this message, as the confidence intervals so obtained will be invalid.")
  
  # -- Compute pos/neg LRs and their BCa confidence intervals -- #
  negLR <- unname( ( 1 - cs[,"sens"] ) / cs[,"spec"]  )
  negLRexact <- unname( ( 1-csExact[,"sens"] ) / csExact[,"spec"] ) # Could also just use csExact[,"negLR"] here, but not in the line above it
  if( all( specb != 0L ) ) {
    negLR.ci <- bca( ( 1 - sensb) / specb, negLR, conf = ci.width, ... )$bca[4:5]
  } else {
    negLR.ci <- 1/bca( specb / ( 1 - sensb), 1/negLR, conf = ci.width, ... )$bca[c(4,5)]
  }
  posLR <- unname( cs[,"sens"] / ( 1 - cs[,"spec"] ) )
  posLRexact <- unname( csExact[,"sens"] / ( 1 - csExact[,"spec"] ) )
  if( all( specb != 1L ) ) {
    posLR.ci <- bca( sensb / ( 1 - specb ), posLR, conf = ci.width, ... )$bca[4:5]
  } else {
    posLR.ci <- 1/bca( ( 1 - specb ) / sensb, 1/posLR, conf = ci.width, ... )$bca[c(5,4)] # Reversed because the order inverts when you take the reciprocal
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
  ci.width = ci.width
  )
}

#' Internal function to draw a set of sensitivities or specificities
#' This is intended for the case where testPos == totalDzPos or testNeg == totalDzNeg.
#' @param n The total number of positives/negatives in the population.
#' @param R is the number of replications in each round of the bootstrap (has been tested at 10,000 or greater).
#' @param consistentQuantile Defaults to 0.5, i.e. the median. Finds the lowest probability for which the random draws are likely to be consistently one, where consistently is defined by this value (i.e. at .5, a simple majority of the time is enough for consistency)
#' @param verbose Whether to display internal operations as they happen.
#' @param parameters List of control parameters (shrink, tol, nEach) for sequential search.
#' @param method Either "deterministic" or "search". The former is faster and more accurate. Thanks to an anonymous reviewer for pointing out the utility of the binomial distribution in solving this problem.
drawMaxedOut <- function( n, R, consistentQuantile = 0.5, verbose, parameters=list(shrink=5,tol=.0005,nEach=80), method = "deterministic" ) {
  # Type checking
  if( class(consistentQuantile) != "numeric" | consistentQuantile > 1 | consistentQuantile < 0 )  stop("consistentQuantile must be a decimal number between 0 and 1")
  # Find lowest probability that consistently produces 1's, where "consistently" is defined by consistentQuantile
  if( method == "search" ) {
    lprb <- sequentialGridSearch( # lowest probability that consistently produces 1's 
      f=identity, # We just want to minimize pr
      constraint=function(probs,...) vapply( probs, FUN=medianConsistentlyOne, FUN.VALUE=NA, consistentQuantile = consistentQuantile, ... ),
      bounds=c(0,1), 
      verbose=verbose,
      size=n, R=R, warn=FALSE,
      shrink=parameters$shrink,
      tol=parameters$tol,
      nEach=parameters$nEach
    )
  }
  else if( method == "deterministic" ) {
    lprb <- exp( log( consistentQuantile ) / n )
  }
  else {
    stop("Method",method,"not supported")
  }
  # Draw samples based on this probability and return result
  res <- stats::rbinom(R, size=n, prob=lprb)/n
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
#' @param digits Number of digits to round to for display purposes
#' @param \dots Pass-alongs (currently ignored).
#' @return Returns x unaltered.
#' @method print lrtest 
#' @export
#' @examples
#' \dontrun{
#' print.lrtest( BayesianLR.test( 500, 500, 300, 500 ), digits = 4 )
#' }
print.lrtest <- function( x, digits = 3, ... ) {
  cat("\n")
  cat("Likelihood ratio test of a 2x2 table")
  cat("\n\n")
  cat("data:\n")
  print(x$inputs)
  cat( sprintf( paste0( "Positive LR: %01.", digits, "f (%01.", digits, "f - %01.", digits, "f)\n" ), x$posLR, x$posLR.ci[1], x$posLR.ci[2] ) )
  cat( sprintf( paste0( "Negative LR: %01.", digits, "f (%01.", digits, "f - %01.", digits, "f)\n" ), x$negLR, x$negLR.ci[1], x$negLR.ci[2] ) )
  cat( paste0( attr(x,"ci.width")*100, "% confidence intervals computed via ", attr(x,"ci.type"), " bootstrapping.\n" ) )
  cat( "Note: This procedure depends on repeated random sampling.  As such it is subject to some variability in results.\n  Variability is minimized by large numbers of replications (typically 10,000), and averaging 5 repeated results,\n but with small sample sizes or sensitivity or specificity near 0 or 1, variability becomes more pronounced.\n  This is not an error, it is a function of the nature of the procedure." )
  invisible(x)
}