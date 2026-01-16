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

#### For testing 
library(readxl); 
if( Sys.info()["user"]=="Ellner") setwd("c:/repos/luckieR/stagecoach"); 
M = as.matrix(read_xlsx("Caswell_P.xlsx")); 
F = as.matrix(read_xlsx("Caswell_B.xlsx")); 

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
        LE = colSums(N); 
  }
  return(LE)
}  

########################################################################
##  Function to calculate modified transition matrix D_i for stage i, ##
##      where individuals entering stage i stay there for one time    ##
##      step and then die. This is a utility for meantime() and       ##
##      generally not of any value to users. 
##  Equation 8 in Cochran and Ellner 1992                             ##
##  Argument/s:                                                       ##
##     M: state transition matrix M                                   ##
##  Return Value/s:                                                   ##
##     Array of D of transition matrices, such that D[,,i] is D_i     ##
########################################################################
Di_mat <- function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  results <- array(0, dim = c(dim(M),ncol(M)))
  for (i  in 1:ncol(M)) {
      results[,,i] = M; 
      results[,i,i] = 0; 
  }
  return(results) 
}

############################################################################
##  Function to compute the average time to reach stage i from stage j,   ##
##         conditional on that happening                                  ## 
##  Equation 9 in Cochran and Ellner 1992                                 ##
##  Argument/s:                                                           ##
##     M: state transition matrix M                                       ##
##  Return Value/s:                                                       ##
##     Matrix, (i,j) entry = conditional mean time to reach i from j      ##
############################################################################

meantime <- function(M){
  if(!is.matrix(M)) {stop("This is not a matrix")}
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  D <- Di_mat(M)
  
  results <- matrix(0, nrow(M), ncol= ncol(M))
  # for loop for each Di_mat
  for (i in 1:dim(results)[1]) {
    Imat <- diag(x=1, dim(D[,,i])[1])
    den <- solve(Imat - D[,,i])  
    num <- den %*% den
    results[i,] <- (num[i,]/den[i,])
  }

  results[is.nan(results)] <- NA
  # NA used because can't get to this state from previous state
  
  return(results)
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
  if(nrow(M)!=ncol(M)) {stop("This is not a square matrix")}
  if(!is.numeric(M)) {stop("This is not numeric")}
  
  LE <- life_expectancy(M)
 
 # must add each stage number to the correct vector
  MT <- meantime(M)
  results <- matrix(0,nrow = nrow(M),ncol = ncol(M))
  for (i in 1:nrow(M)) {
    for (j in 1:nrow(M)) {
      results[i,j]  = MT[i,j] + LE[i] + 1
      if(i < j){
         results[i,j] = -1
       } else{
         results[i,j] = MT[i,j] + LE[i] + 1
       }
    }
  }
  results[results == 0] = NA
  return(results)
}

