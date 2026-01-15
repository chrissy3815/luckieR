###########################################################################
## R translation of formulas for mean/sd of conditional lifespan,from
##     Cochran, M. E. and S. Ellner.  1992.  Simple methods for 
##     calculating age-based life history parameters for stage-
##     structured populations. Ecological Monographs 62:345-364.
##
## Original coding by Erin Souder, merge into luckieR by Steve Ellner 
###########################################################################

######################### NOTATION CONVENTIONS ############################
##      luckieR M = StageCoach P (state transitions) 
##      luckieR F = StageCoach B (sexual births) + F (clonal births, fission) 
########################################################################### 


########################################################################### 
## Compute mean life expectancy given current state, for all states,     ##
##		given by row sums of the fundamental matrix                      ##
##  Arguments:                                                           ##
##  		M: state transition matrix M                                 ##
##  Return Value/s:                                                      ##
##        Mean lifespan expectancy for each state i  					 ##
########################################################################### 
life_expectancy = function(M){ 
  matDim = nrow(M); 
  N <- try(solve(diag(matDim) - M), silent = TRUE)
  if (inherits(N, "try-error")) {
        LE = rep(NA,matDim)
		} else {
        LE <- rep(1, matDim) %*% N
  }
  return(LE)
}  


#########################################################################
##  Function to determine mean of lifespan conditional on individual   ##
##      reaching state i before death.                                 ##
##  Equation 6 in Cochran and Ellner 1992                              ##
##  Arguments:                                                         ##
##  		M: state transition matrix M                               ##
##  Return Value/s:                                                    ##
##        Conditional mean lifespan for each state i                   ##
#########################################################################
conditional_mean_lifeSpan <- function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(!is.square.matrix(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  LE <- life_expectancy(M)
 
 # must add each stage number to the correct vector
  MT <- meantime(M)
  results <- matrix(0,nrow = nrow(M),ncol = ncol(M))
  for (i in 1:nrow(M)) {
    for (j in 1:nrow(M) {
      results[i,j]  = MT[i,j] + LE[i] + 1
      if(i < j){
         results[i,j] = -1
       } else{
         results[i,j] = MT[i,j] + LE[i] + 1
       }
    }
  }
  results[results == 0] <- NA
  return(results)
}

