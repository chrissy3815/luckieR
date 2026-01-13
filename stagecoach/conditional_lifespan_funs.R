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

conditional_mean_lifeSpan <- function(M){
  
#########################################################################
##  Function to determine mean of lifespan conditional on individual   ##
##      reaching state i before death.                                 ##
##  Equation 6 in Cochran and Ellner 1992                              ##
##  Arguments:                                                         ##
##  		M: state transition matrix M                               ##
##  Return Value/s:                                                    ##
##        Conditional mean lifespan for each state i                   ##
#########################################################################

  P = M; ## go to StageCoach notation  
  if(!is.matrix(P)) {stop("This is not a matrix")}
  if(!is.square.matrix(P)) {stop("This is not a square matrix")}
  if(!is.numeric(P)) {stop("This is not numeric")}
  

  ## Compute life expectancy given current state 
  matDim = nrow(P); 
  N <- try(solve(diag(matDim) - P), silent = TRUE)
  if (inherits(N, "try-error")) {
        LE = rep(NA,matDim)
		} else {
        LE <- rep(1, matDim) %*% N
    }

################ STEVE STOPPED HERE   
  
  LE <- life_expectancy(P)
 # must add each stage number to the correct vector
  MT <- meantime(P)
  results <- matrix(0,nrow = nrow(P),ncol = ncol(P))
  for (i in 1:dim(P)[1]) {
   for (j in 1:dim(P)[2]) {
     results[i,j] <- MT[i,j] + LE[i] + 1
     if(i < j){
       results[i,j] <- -1
     }
     else{
       results[i,j] <- MT[i,j] + LE[i] + 1
   }
     }
 }
  results[results < 0] < NA
  return(results)
}